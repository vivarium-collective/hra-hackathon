"""Run all HRA hackathon models through COPASI via process-bigraph."""

from __future__ import annotations

import json
import os
from pathlib import Path

from process_bigraph import allocate_core, Composite, gather_emitter_results
from process_bigraph.emitter import emitter_from_wires

from processes.copasi_process import CopasiUTCStep

from experiments.biomodels_fetch import get_metadata_utf8, load_biomodel
from experiments.model_ids import HRA_MODEL_IDS
from experiments.species_names import resolve_species_names
from experiments.sbml_utils import get_sbml_info
from experiments.plotting import extract_result, make_plotly_figure, build_results_page

OUTDIR = "out_hra_models"


def make_document(model_id, sbml_path, duration, n_points):
    """Build a process-bigraph document with a COPASI step and RAM emitter."""
    step_key = f"{model_id}_copasi"
    return {
        "schema": {
            "species_concentrations": "map[float]",
            "species_counts": "map[float]",
        },
        "state": {
            "species_concentrations": {},
            "species_counts": {},
            "results": {},
            f"{step_key}_step": {
                "_type": "step",
                "address": "local:CopasiUTCStep",
                "config": {
                    "model_source": sbml_path,
                    "time": duration,
                    "n_points": n_points,
                },
                "inputs": {
                    "species_concentrations": ["species_concentrations"],
                    "species_counts": ["species_counts"],
                },
                "outputs": {
                    "result": ["results", step_key],
                },
            },
            "emitter": emitter_from_wires({
                "results": ["results"],
                "species_concentrations": ["species_concentrations"],
            }),
        },
    }


def run_document(document, core, name):
    """Run a composite document, return emitter data."""
    os.makedirs(OUTDIR, exist_ok=True)

    duration = max(
        (node.get("config", {}).get("time", 0)
         for node in document["state"].values()
         if isinstance(node, dict) and node.get("_type") == "step"),
        default=10.0,
    )

    sim = Composite(document, core=core)
    print(f"  Running {name} for {duration}s ...")
    sim.run(duration)
    print(f"  Done: {name}")

    Path(os.path.join(OUTDIR, f"{name}.json")).write_text(
        json.dumps(document, indent=2), encoding="utf-8")
    try:
        serialized = core.serialize(sim.schema, sim.state)
        Path(os.path.join(OUTDIR, f"{name}_state.json")).write_text(
            json.dumps(serialized, indent=2), encoding="utf-8")
    except Exception as e:
        print(f"  WARNING: Could not serialize state: {e}")

    return gather_emitter_results(sim).get(("emitter",), [])


def run_hra_models(core, model_ids=None):
    """Run simulations for the given BioModels IDs (defaults to HRA_MODEL_IDS)."""
    if model_ids is None:
        model_ids = HRA_MODEL_IDS

    core.register_link("CopasiUTCStep", CopasiUTCStep)
    core.register_type("numeric_result", {
        "time": "list[float]",
        "columns": "list[string]",
        "values": "list[list[float]]",
    })

    default_duration, default_n_points = 10.0, 100
    status = []
    model_results = []

    for model_id in model_ids:
        print(f"\n{'='*60}\nProcessing {model_id}\n{'='*60}")
        try:
            metadata = get_metadata_utf8(model_id)
            try:
                sbml_path, duration, n_points = load_biomodel(model_id, metadata)
            except ValueError:
                sbml_path = os.path.abspath(f"models/{model_id}.xml")
                if not os.path.exists(sbml_path):
                    raise
                print(f"  No SED-ML found, using local SBML with defaults")
                duration, n_points = default_duration, default_n_points

            print(f"  SBML: {sbml_path}  Duration: {duration}s  Points: {n_points}")

            doc = make_document(model_id, sbml_path, duration, n_points)
            emitter_data = run_document(doc, core, f"{model_id}_copasi")

            step_key = f"{model_id}_copasi"
            result_data = extract_result(emitter_data, step_key)
            info = get_sbml_info(sbml_path)

            if result_data:
                time_unit = info.get("time_unit", "s")
                title = f"{info.get('name', model_id)} — Time Course Simulation"
                print(f"  Resolving species names from SBML annotations...")
                name_map = resolve_species_names(sbml_path)
                fig = make_plotly_figure(result_data, title, time_unit=time_unit,
                                        name_map=name_map)
                plot_html = fig.to_html(full_html=False, include_plotlyjs=False)
                model_results.append({
                    "model_id": model_id,
                    "info": info,
                    "duration": duration,
                    "n_points": n_points,
                    "time_unit": time_unit,
                    "plot_html": plot_html,
                })
            else:
                print(f"  No result data for {step_key}")

            status.append((model_id, True, ""))

        except Exception as e:
            print(f"  ERROR: {e}")
            status.append((model_id, False, str(e)))

    build_results_page(model_results, OUTDIR)

    ok = sum(1 for _, s, _ in status if s)
    print(f"\n{'='*60}\nSUMMARY: {ok}/{len(status)} succeeded\n{'='*60}")
    for mid, success, err in status:
        print(f"  {mid}: {'OK' if success else f'FAIL: {err}'}")


if __name__ == "__main__":
    import sys
    core = allocate_core()
    ids = sys.argv[1:] or None
    run_hra_models(core, model_ids=ids)
