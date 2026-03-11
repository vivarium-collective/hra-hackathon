"""Run all HRA hackathon models through COPASI via process-bigraph."""

from __future__ import annotations

import builtins
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Any

import time as _time

import biomodels
import libsbml
import libsedml
import plotly.graph_objects as go
import requests
from process_bigraph import allocate_core, Composite, gather_emitter_results
from process_bigraph.emitter import emitter_from_wires

from processes.copasi_process import CopasiUTCStep

HRA_MODEL_IDS = [
    "BIOMD0000000356",
    "BIOMD0000000620",
    "BIOMD0000000621",
    "BIOMD0000000833",
    "BIOMD0000000854",
    "MODEL1209260000",
    "MODEL1912090001",
    "MODEL2401110001",
]

OUTDIR = "out_hra_models"

_SBML_RE = re.compile(r"\.(xml|sbml)$", re.IGNORECASE)
_SEDML_RE = re.compile(r"\.sedml$", re.IGNORECASE)


# --- BioModels file helpers ---

def _file_name(obj):
    return getattr(obj, "name", str(obj))


def _iter_entry_files(entry):
    if isinstance(entry, (list, tuple)):
        return entry
    if isinstance(entry, dict):
        for key in ("files", "main_files", "model_files"):
            v = entry.get(key)
            if isinstance(v, (list, tuple)):
                return v
        return []
    try:
        return list(entry)
    except TypeError:
        return []


def _find_file(entry_files, regex, prefer_keywords=None):
    candidates = [f for f in entry_files if regex.search(_file_name(f))]
    if prefer_keywords:
        for kw in prefer_keywords:
            for c in candidates:
                if kw in _file_name(c).lower():
                    return c
    return candidates[0] if candidates else None


def _fetch_file(file_entry, out_dir):
    f = biomodels.get_file(file_entry)
    if isinstance(f, (str, os.PathLike)) and os.path.exists(str(f)):
        return str(f)
    out_path = os.path.join(out_dir, _file_name(file_entry))
    if isinstance(f, bytes):
        Path(out_path).write_bytes(f)
    else:
        Path(out_path).write_text(str(f), encoding="utf-8")
    return out_path


def _get_metadata_utf8(model_id):
    """biomodels.get_metadata wrapper — works around ASCII encoding bug in cache."""
    try:
        return biomodels.get_metadata(model_id)
    except UnicodeDecodeError:
        orig = builtins.open
        def patched(path, *a, **kw):
            if str(path).endswith(model_id):
                kw.setdefault("encoding", "utf-8")
            return orig(path, *a, **kw)
        builtins.open = patched
        try:
            return biomodels.get_metadata(model_id)
        finally:
            builtins.open = orig


# --- SED-ML parsing ---

def _extract_utc(sedml_path):
    """Parse SED-ML and return (duration, n_points) from first UniformTimeCourse."""
    doc = libsedml.readSedMLFromFile(str(sedml_path))
    if doc is None or doc.getNumErrors() > 0:
        raise RuntimeError(f"SED-ML parse error: {sedml_path}")
    for i in range(doc.getNumSimulations()):
        sim = doc.getSimulation(i)
        if sim and all(hasattr(sim, m) for m in
                       ("getInitialTime", "getOutputStartTime", "getOutputEndTime", "getNumberOfPoints")):
            duration = float(sim.getOutputEndTime()) - float(sim.getOutputStartTime())
            return duration, int(sim.getNumberOfPoints())
    raise ValueError("No UniformTimeCourse found in SED-ML")


def _resolve_sbml_from_sedml(sedml_path, fallback_sbml):
    """If SED-ML references a local SBML file, resolve it; otherwise use fallback."""
    doc = libsedml.readSedMLFromFile(str(sedml_path))
    if doc and doc.getNumModels() > 0:
        src = doc.getModel(0).getSource() or ""
        if src and not src.startswith(("http", "urn:", "biomodels:", "BIOMD")):
            candidate = os.path.join(os.path.dirname(os.path.abspath(sedml_path)), src)
            if os.path.exists(candidate):
                return candidate
    return fallback_sbml


# --- Load a BioModel ---

def load_biomodel(model_id, metadata):
    """Fetch SED-ML + SBML from BioModels, return (sbml_path, duration, n_points)."""
    files = list(_iter_entry_files(metadata))
    sedml_files = [f for f in files if _SEDML_RE.search(_file_name(f))]
    sbml_files = [f for f in files
                  if _SBML_RE.search(_file_name(f)) and not _SEDML_RE.search(_file_name(f))]

    if not sedml_files or not sbml_files:
        raise ValueError(f"{model_id}: missing SED-ML or SBML in BioModels entry")

    with tempfile.TemporaryDirectory(prefix=f"biomodel_{model_id}_") as tmp:
        sedml_path = _fetch_file(sedml_files[0], tmp)
        sbml_path = _fetch_file(sbml_files[0], tmp)
        duration, n_points = _extract_utc(sedml_path)
        sbml_path = _resolve_sbml_from_sedml(sedml_path, sbml_path)

        # Copy to stable location
        stable_dir = os.path.join("models", model_id)
        os.makedirs(stable_dir, exist_ok=True)
        stable_sbml = os.path.join(stable_dir, os.path.basename(sbml_path))
        Path(stable_sbml).write_bytes(Path(sbml_path).read_bytes())

    return stable_sbml, duration, n_points


# --- Document creation & running ---

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

    # Infer duration from step config
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

    # Save document and state
    Path(os.path.join(OUTDIR, f"{name}.json")).write_text(
        json.dumps(document, indent=2), encoding="utf-8")
    try:
        serialized = core.serialize(sim.schema, sim.state)
        Path(os.path.join(OUTDIR, f"{name}_state.json")).write_text(
            json.dumps(serialized, indent=2), encoding="utf-8")
    except Exception as e:
        print(f"  WARNING: Could not serialize state: {e}")

    return gather_emitter_results(sim).get(("emitter",), [])


# --- Species name resolution ---

_IDENT_RE = re.compile(r'identifiers\.org/(uniprot|CHEBI|chebi|kegg\.compound)/([A-Za-z0-9:_.-]+)')

# Global cache so we don't re-fetch the same protein across models
_NAME_CACHE: dict[str, str] = {}


def _extract_species_identifiers(sbml_path):
    """Parse SBML annotations and return {species_name: {db: [ids]}}."""
    doc = libsbml.readSBML(sbml_path)
    model = doc.getModel() if doc else None
    if not model:
        return {}

    result = {}
    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        name = sp.getName() or sp.getId()
        annotation = sp.getAnnotationString() or ""
        ids: dict[str, list[str]] = {}
        for match in _IDENT_RE.finditer(annotation):
            db = match.group(1).lower()
            acc = match.group(2)
            ids.setdefault(db, []).append(acc)
        if ids:
            result[name] = ids
    return result


def _fetch_uniprot_name(accession):
    """Fetch recommended protein name from UniProt REST API."""
    if accession in _NAME_CACHE:
        return _NAME_CACHE[accession]
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{accession}.json",
            timeout=10,
        )
        if resp.status_code == 200:
            data = resp.json()
            desc = data.get("proteinDescription", {})
            rec = desc.get("recommendedName", {})
            # Prefer short name, fall back to full name
            short = rec.get("shortNames", [{}])
            if short and short[0].get("value"):
                name = short[0]["value"]
            else:
                name = rec.get("fullName", {}).get("value", "")
            if name:
                _NAME_CACHE[accession] = name
                return name
    except requests.RequestException:
        pass
    _NAME_CACHE[accession] = ""
    return ""


def _fetch_kegg_name(compound_id):
    """Fetch compound name from KEGG REST API."""
    key = f"kegg:{compound_id}"
    if key in _NAME_CACHE:
        return _NAME_CACHE[key]
    try:
        resp = requests.get(f"https://rest.kegg.jp/get/{compound_id}", timeout=10)
        if resp.status_code == 200:
            for line in resp.text.splitlines():
                if line.startswith("NAME"):
                    name = line.split(None, 1)[1].rstrip(";").strip()
                    _NAME_CACHE[key] = name
                    return name
    except requests.RequestException:
        pass
    _NAME_CACHE[key] = ""
    return ""


def _clean_species_name(name):
    """Clean up a species name: strip x##_ prefixes, replace underscores."""
    # Strip leading x##_ prefix (e.g. x9_IRS_1 -> IRS 1)
    cleaned = re.sub(r"^x\d+_", "", name)
    # Strip BioNetGen-style state annotations: "IR(NPXY,Y999~u,...)" -> "IR"
    cleaned = re.sub(r"\(.*\)$", "", cleaned)
    # Replace underscores with spaces
    cleaned = cleaned.replace("_", " ")
    return cleaned.strip() or name


def _resolve_species_names(sbml_path):
    """Build a display name mapping for species in an SBML model.

    Returns {copasi_display_name: human_readable_name}.
    """
    doc = libsbml.readSBML(sbml_path)
    model = doc.getModel() if doc else None
    if not model:
        return {}

    identifiers = _extract_species_identifiers(sbml_path)
    name_map = {}
    uniprot_batch = {}  # name -> accession, for batch lookup

    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        sp_name = sp.getName() or sp.getId()
        sp_ids = identifiers.get(sp_name, {})

        # Collect UniProt IDs for batch fetching
        if "uniprot" in sp_ids:
            acc = sp_ids["uniprot"][0]
            if acc not in _NAME_CACHE:
                uniprot_batch[sp_name] = acc

    # Batch fetch UniProt names (with rate limiting)
    for sp_name, acc in uniprot_batch.items():
        _fetch_uniprot_name(acc)
        _time.sleep(0.1)  # be polite to UniProt API

    # Now build the final name map
    seen_names: dict[str, int] = {}  # track duplicates
    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        sp_name = sp.getName() or sp.getId()
        sp_ids = identifiers.get(sp_name, {})

        resolved = ""
        db_label = ""

        # Try UniProt first
        if "uniprot" in sp_ids:
            acc = sp_ids["uniprot"][0]
            resolved = _NAME_CACHE.get(acc, "")
            if resolved:
                db_label = f"UniProt:{acc}"

        # Try KEGG
        if not resolved and "kegg.compound" in sp_ids:
            kid = sp_ids["kegg.compound"][0]
            resolved = _fetch_kegg_name(kid)
            if resolved:
                db_label = f"KEGG:{kid}"
            _time.sleep(0.1)

        # Use SBML name attribute if it's informative (not same as ID)
        if not resolved:
            sbml_name = sp.getName() or ""
            sbml_id = sp.getId() or ""
            if sbml_name and sbml_name != sbml_id:
                resolved = _clean_species_name(sbml_name)
            else:
                resolved = _clean_species_name(sbml_id)

        # Disambiguate duplicate display names (e.g. multiple "INSR" states)
        count = seen_names.get(resolved, 0)
        seen_names[resolved] = count + 1
        if count > 0:
            resolved = f"{resolved} ({sp_name})"

        name_map[sp_name] = {"display": resolved, "db_id": db_label}

    return name_map


# --- Model info & plotting ---

# Time units known from publications but missing from SBML metadata
_TIME_UNIT_OVERRIDES = {
    "Nyman2011":  "min",   # Dalla Man glucose model family uses minutes
    "Koenig2012": "min",   # hepatic glucose metabolism, rates in µmol/kg/min
}


def _sbml_time_unit(model):
    """Determine the human-readable time unit from an SBML model."""
    tu_id = model.getTimeUnits() or ""
    # Check built-in unit names first
    if tu_id in ("second", "s"):
        return "s"
    if tu_id in ("minute", "min"):
        return "min"
    if tu_id in ("hour", "h"):
        return "h"
    if tu_id in ("day", "d"):
        return "day"

    # Look up the unit definition and compute total multiplier in seconds
    ud = model.getUnitDefinition(tu_id) if tu_id else None
    if ud and ud.getNumUnits() > 0:
        total_seconds = 1.0
        for i in range(ud.getNumUnits()):
            u = ud.getUnit(i)
            if u.getKind() == libsbml.UNIT_KIND_SECOND:
                total_seconds = u.getMultiplier() * (10 ** u.getScale()) ** u.getExponent()
                break
        if abs(total_seconds - 1.0) < 0.1:
            return "s"
        if abs(total_seconds - 60.0) < 1.0:
            return "min"
        if abs(total_seconds - 3600.0) < 60.0:
            return "h"
        if abs(total_seconds - 86400.0) < 600.0:
            return "day"

    # Fallback: check substance/time unit definitions for clues
    for i in range(model.getNumUnitDefinitions()):
        ud = model.getUnitDefinition(i)
        uid = ud.getId().lower()
        if "time" in uid or uid in ("minute", "min", "second", "day", "hour"):
            for j in range(ud.getNumUnits()):
                u = ud.getUnit(j)
                if u.getKind() == libsbml.UNIT_KIND_SECOND:
                    secs = u.getMultiplier() * (10 ** u.getScale()) ** u.getExponent()
                    if abs(secs - 60.0) < 1.0:
                        return "min"
                    if abs(secs - 86400.0) < 600.0:
                        return "day"
                    if abs(secs - 3600.0) < 60.0:
                        return "h"
    # Check manual overrides by model name
    model_name = (model.getName() or model.getId() or "")
    for key, unit in _TIME_UNIT_OVERRIDES.items():
        if key.lower() in model_name.lower():
            return unit

    return "s"  # default


def _get_sbml_info(sbml_path):
    """Extract model metadata from SBML file."""
    doc = libsbml.readSBML(sbml_path)
    m = doc.getModel() if doc else None
    if not m:
        return {}
    compartments = [m.getCompartment(i).getId() for i in range(m.getNumCompartments())]
    return {
        "name": (m.getName() or m.getId() or "").replace("_", " "),
        "num_species": m.getNumSpecies(),
        "num_reactions": m.getNumReactions(),
        "num_compartments": m.getNumCompartments(),
        "compartments": compartments,
        "num_parameters": m.getNumParameters(),
        "time_unit": _sbml_time_unit(m),
    }


def _extract_result(emitter_data, step_key):
    """Pull the time-course result dict from emitter snapshots."""
    for snapshot in emitter_data:
        if snapshot and step_key in snapshot.get("results", {}):
            r = snapshot["results"][step_key]
            if r:
                return r
    return None


def _downsample(time_points, values, max_points=100):
    """Downsample time series to at most max_points, keeping first and last."""
    n = len(time_points)
    if n <= max_points:
        return time_points, values
    step = max(1, n // max_points)
    indices = list(range(0, n, step))
    if indices[-1] != n - 1:
        indices.append(n - 1)
    return [time_points[i] for i in indices], [values[i] for i in indices]


def _make_plotly_figure(result_data, title, time_unit="s", name_map=None):
    """Build a Plotly figure from time-course result data."""
    time_pts, vals = _downsample(result_data["time"], result_data["values"])
    # Round to 6 significant figures to reduce HTML size
    time_pts = [round(t, 6) for t in time_pts]
    vals = [[float(f"{v:.6g}") for v in row] for row in vals]

    fig = go.Figure()
    for i, species in enumerate(result_data["columns"]):
        entry = (name_map or {}).get(species, {})
        display = entry.get("display", species) if entry else species
        db_id = entry.get("db_id", "") if entry else ""
        db_line = f"<br>{db_id}" if db_id else ""
        fig.add_trace(go.Scatter(
            x=time_pts,
            y=[row[i] for row in vals],
            mode="lines",
            name=display,
            hovertemplate=f"<b>{display}</b>{db_line}<br>"
                          f"Time: %{{x:.2f}}<br>Conc: %{{y:.4g}}<extra></extra>",
        ))
    unit_labels = {"s": "seconds", "min": "minutes", "h": "hours", "day": "days"}
    fig.update_layout(
        title=title,
        xaxis_title=f"Time ({unit_labels.get(time_unit, time_unit)})",
        yaxis_title="Concentration",
        hovermode="x unified",
        legend=dict(title="Species", font=dict(size=10)),
        template="plotly_white",
        height=500,
    )
    return fig


def build_results_page(model_results, outdir):
    """Generate a single HTML page with all model results."""
    sections = []
    for entry in model_results:
        info = entry["info"]
        model_id = entry["model_id"]
        name = info.get("name", model_id)
        biomodels_url = f"https://www.ebi.ac.uk/biomodels/{model_id}"

        # Info card
        sections.append(f"""
        <div class="model-card" id="{model_id}">
            <h2>{name}</h2>
            <table class="info-table">
                <tr><td>BioModels ID</td><td><a href="{biomodels_url}" target="_blank">{model_id}</a></td></tr>
                <tr><td>Species</td><td>{info.get('num_species', '?')}</td></tr>
                <tr><td>Reactions</td><td>{info.get('num_reactions', '?')}</td></tr>
                <tr><td>Compartments</td><td>{', '.join(info.get('compartments', [])) or '?'}</td></tr>
                <tr><td>Parameters</td><td>{info.get('num_parameters', '?')}</td></tr>
                <tr><td>Simulation</td><td>{entry['duration']} {entry.get('time_unit', 's')}, {entry['n_points']} time points</td></tr>
            </table>
            <div class="plot">{entry['plot_html']}</div>
        </div>
        """)

    # Table of contents
    toc_items = "".join(
        f'<li><a href="#{e["model_id"]}">{e["info"].get("name", e["model_id"])}</a></li>'
        for e in model_results
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>HRA Hackathon — Simulation Results</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         max-width: 1100px; margin: 0 auto; padding: 20px; background: #fafafa; color: #333; }}
  h1 {{ border-bottom: 2px solid #2c5f8a; padding-bottom: 10px; color: #2c5f8a; }}
  h2 {{ color: #2c5f8a; margin-top: 0; }}
  nav ul {{ columns: 2; list-style: none; padding: 0; }}
  nav li {{ padding: 4px 0; }}
  nav a {{ color: #2c5f8a; text-decoration: none; }}
  nav a:hover {{ text-decoration: underline; }}
  .model-card {{ background: white; border-radius: 8px; padding: 24px;
                 margin: 24px 0; box-shadow: 0 1px 4px rgba(0,0,0,0.1); }}
  .info-table {{ border-collapse: collapse; margin: 12px 0 20px; }}
  .info-table td {{ padding: 4px 16px 4px 0; }}
  .info-table td:first-child {{ font-weight: 600; color: #555; white-space: nowrap; }}
  .info-table a {{ color: #2c5f8a; }}
  .plot {{ margin-top: 12px; }}
</style>
</head>
<body>
<h1>HRA Hackathon — Simulation Results</h1>
<p>{len(model_results)} models simulated with <a href="https://copasi.org/">COPASI</a>
via <a href="https://github.com/vivarium-collective/process-bigraph">process-bigraph</a>.
Hover over traces for values, click legend entries to toggle species, drag to zoom.</p>
<nav><ul>{toc_items}</ul></nav>
{"".join(sections)}
</body>
</html>"""

    filepath = os.path.join(outdir, "results.html")
    os.makedirs(outdir, exist_ok=True)
    Path(filepath).write_text(html, encoding="utf-8")
    print(f"\nResults page saved: {filepath}")
    return filepath


# --- Main ---

def run_hra_models(core):
    core.register_link("CopasiUTCStep", CopasiUTCStep)
    core.register_type("numeric_result", {
        "time": "list[float]",
        "columns": "list[string]",
        "values": "list[list[float]]",
    })

    default_duration, default_n_points = 10.0, 100
    status = []
    model_results = []

    for model_id in HRA_MODEL_IDS:
        print(f"\n{'='*60}\nProcessing {model_id}\n{'='*60}")
        try:
            metadata = _get_metadata_utf8(model_id)
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
            result_data = _extract_result(emitter_data, step_key)
            info = _get_sbml_info(sbml_path)

            if result_data:
                time_unit = info.get("time_unit", "s")
                title = f"{info.get('name', model_id)} — Time Course Simulation"
                print(f"  Resolving species names from SBML annotations...")
                name_map = _resolve_species_names(sbml_path)
                fig = _make_plotly_figure(result_data, title, time_unit=time_unit,
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
    core = allocate_core()
    run_hra_models(core)
