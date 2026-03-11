import os
from pathlib import Path
from typing import Dict, Any

from pandas import DataFrame
from process_bigraph import Process, Step, Composite, gather_emitter_results
import COPASI
from basico import (
    load_model,
    get_species,
    get_reactions,
    set_species,
    run_time_course,
    run_steadystate,
)
import COPASI

from processes import model_path_resolution


def _set_initial_concentrations(changes, dm):
    """
    changes: iterable of (species_name, value) pairs
    dm: COPASI DataModel as returned by basico.load_model
    """
    model = dm.getModel()
    assert isinstance(model, COPASI.CModel)

    references = COPASI.ObjectStdVector()

    for name, value in changes:
        species = model.getMetabolite(name)
        if species is None:
            print(f"Species {name} not found in model")
            continue
        assert isinstance(species, COPASI.CMetab)
        species.setInitialConcentration(float(value))
        references.append(species.getInitialConcentrationReference())

    if len(references) > 0:
        model.updateInitialValues(references)


def _get_transient_concentration(name, dm):
    """
    Return the *current* concentration (not initial) of a species.
    """
    model = dm.getModel()
    assert isinstance(model, COPASI.CModel)

    species = model.getMetabolite(name)
    if species is None:
        print(f"Species {name} not found in model")
        return None
    assert isinstance(species, COPASI.CMetab)
    return float(species.getConcentration())

class BaseCopasi:
    cmodel = None
    dm = None
    species_ids = None
    reaction_ids = None
    sbml_to_name = None

    def interpret_sbml(self):
        model_source = self.config['model_source']

        # ---- Load COPASI model ----
        self.dm = load_model(model_path_resolution(model_source))
        if self.dm is None:
            raise RuntimeError(
                f"load_model({model_source!r}) returned None. "
                "Check that the file exists and is a valid COPASI/SBML model."
            )

        self.cmodel = self.dm.getModel()

        spec_df = get_species(model=self.dm)

        # External canonical IDs: SBML IDs
        self.species_ids = spec_df["sbml_id"].tolist()

        # Mapping: SBML ID -> COPASI display name (index)
        self.sbml_to_name = {
            spec_df.loc[name, "sbml_id"]: name
            for name in spec_df.index
        }

        rxn_df = get_reactions(model=self.dm)
        # These are typically SBML reaction ids already
        self.reaction_ids = rxn_df.index.tolist()

    def get_concentrations_from_sbml(self) -> Dict[str, Any]:
        return {
            "species_concentrations": {
                sbml_id: _get_transient_concentration(
                    name=self.sbml_to_name[sbml_id],  # COPASI name
                    dm=self.dm
                )
                for sbml_id in self.species_ids
            }
        }




class CopasiUTCStep(Step, BaseCopasi):

    config_schema = {
        'model_source': 'string',
        'time': 'float',
        'n_points': 'integer',
        'method': 'string',
        'r_tol': 'float',
        'a_tol': 'float',
    }

    def initialize(self, config=None):
        self.interpret_sbml()

        self.interval = float(self.config.get('time', 1.0))
        self.n_points = int(self.config.get('n_points', 2))
        if self.n_points < 2:
            raise ValueError("n_points must be >= 2")
        self.intervals = self.n_points - 1

    def initial_state(self) -> Dict[str, Any]:
        return self.get_concentrations_from_sbml()

    def inputs(self):
        return {
            'species_concentrations': 'map[float]',
            'species_counts': 'map[float]',
        }

    def outputs(self):
        return {
            'result': 'numeric_result',
        }

    def _run_with_roadrunner(self) -> DataFrame:
        """Fallback simulator for models COPASI can't handle."""
        import tellurium as te
        rr = te.loadSBMLModel(model_path_resolution(self.config['model_source']))
        result = rr.simulate(0, self.config['time'], self.n_points)
        cols = [c for c in result.colnames if c != 'time']
        clean_cols = [c.strip('[]').split('].')[-1] for c in cols]
        tc = DataFrame(
            data=result[:, 1:],
            index=result[:, 0],
            columns=clean_cols,
        )
        tc.index.name = 'Time'
        return tc

    def update(self, inputs):
        # Apply incoming concentrations
        spec_data = inputs.get('counts', {}) or {}
        changes = [
            (name, float(value))
            for name, value in spec_data.items()
            if name in self.species_ids
        ]

        if changes:
            _set_initial_concentrations(changes, self.dm)

        # --- Run COPASI time course ---
        tc_kwargs = dict(
            start_time=0.0,
            duration=self.config['time'],
            intervals=self.intervals,
            update_model=True,
            use_sbml_id=True,
            model=self.dm,
        )
        if self.config.get('method'):
            tc_kwargs['method'] = self.config['method']
        if self.config.get('r_tol'):
            tc_kwargs['r_tol'] = self.config['r_tol']
        if self.config.get('a_tol'):
            tc_kwargs['a_tol'] = self.config['a_tol']

        tc: DataFrame = run_time_course(**tc_kwargs)

        # If COPASI fails (NaN in results), fall back to RoadRunner
        if tc.isnull().any().any():
            print(f"  COPASI produced NaN, falling back to RoadRunner")
            tc = self._run_with_roadrunner()

        time_list = tc.index.to_list()

        result = {
            "time": time_list,
            "columns": [c for c in tc.columns],
            "values": tc.values.tolist(),
        }

        return {"result": result}



class CopasiSteadyStateStep(Step, BaseCopasi):

    config_schema = {
        'model_source': 'string',
        'time': 'float',  # kept for symmetry, not used
    }

    def initialize(self, config=None):
        self.interpret_sbml()

    # ------------------------------------------------
    # initial state (SBML IDs externally)
    # ------------------------------------------------
    def initial_state(self) -> Dict[str, Any]:
        return self.get_concentrations_from_sbml()

    # ------------------------------------------------
    # ports
    # ------------------------------------------------
    def inputs(self):
        # Externally everything uses SBML IDs
        return {
            'species_concentrations': 'map[float]',  # SBML IDs
            'counts': 'map[float]',          # SBML IDs
        }

    def outputs(self):
        # Match TelluriumSteadyStateStep: nested results
        return {
            'results': 'any',
        }

    # ------------------------------------------------
    # steady-state update
    # ------------------------------------------------
    def update(self, inputs):
        # 1) Prefer counts, otherwise concentrations (keys are SBML IDs)
        spec_data = (
            inputs.get('counts')
            or inputs.get('concentrations')
            or {}
        )

        # Convert SBML IDs -> COPASI names for internal set
        changes = []
        for sbml_id, value in spec_data.items():
            name = self.sbml_to_name.get(sbml_id)
            if name is not None:
                changes.append((name, float(value)))

        if changes:
            _set_initial_concentrations(changes, self.dm)

        # 2) Run COPASI steady-state task
        # (use_sbml_id affects task I/O naming, but we read from get_species anyway)
        run_steadystate(
            update_model=True,
            use_sbml_id=True,
            model=self.dm,
        )

        # 3) Read back steady-state species concentrations (SBML IDs externally)
        spec_df = get_species(model=self.dm)
        # spec_df is indexed by COPASI name, with 'sbml_id' and 'concentration' columns
        species_conc_ss = {}
        for name in spec_df.index:
            sbml_id = spec_df.loc[name, "sbml_id"]
            if sbml_id in self.species_ids:
                species_conc_ss[sbml_id] = float(spec_df.loc[name, "concentration"])

        # 4) Steady-state reaction fluxes
        rxn_df = get_reactions(model=self.dm)
        reaction_fluxes_ss = {
            rid: float(rxn_df.loc[rid, 'flux'])
            for rid in self.reaction_ids
            if rid in rxn_df.index
        }

        # 5) Package as one-point "time series" (t = 0.0) to match Tellurium
        time_list = [0.0]

        species_json = {sid: [val] for sid, val in species_conc_ss.items()}
        flux_json = {rid: [val] for rid, val in reaction_fluxes_ss.items()}

        results = {
            "time": time_list,
            "species_concentrations": species_json,  # SBML IDs as keys
            "fluxes": flux_json,
        }

        return {"results": results}


class CopasiUTCProcess(Process, BaseCopasi):

    config_schema = {
        'model_source': 'string',
        'time': 'float',
        'intervals': 'integer',
    }

    def initialize(self, config=None):
        self.interpret_sbml()

        # ---- Sim parameters ----
        self.time = float(self.config.get("time", 1.0))
        self.intervals = int(self.config.get("intervals", 10))

    # -----------------------------------------------------------------
    # initial state
    # -----------------------------------------------------------------
    def initial_state(self) -> Dict[str, Any]:
        return self.get_concentrations_from_sbml()

    # -----------------------------------------------------------------
    # I/O schema
    # -----------------------------------------------------------------
    def inputs(self):
        return {
            "species_concentrations": "map[float]",  # SBML IDs
            "species_counts": "map[float]",          # SBML IDs
        }

    def outputs(self):
        return {
            "species_concentrations": "map[float]",  # SBML IDs
            "fluxes": "map[float]",
            "time": "list[float]",
        }

    # -----------------------------------------------------------------
    # update
    # -----------------------------------------------------------------
    def update(self, inputs, interval):
        # --- 1) Determine incoming species map (SBML IDs)
        incoming = (
            inputs.get("species_counts")
            or inputs.get("species_concentrations")
            or {}
        )

        # Convert SBML IDs → COPASI names for internal setting
        changes = []
        for sbml_id, value in incoming.items():
            name = self.sbml_to_name.get(sbml_id)
            if name is not None:
                changes.append((name, float(value)))

        if changes:
            _set_initial_concentrations(changes, self.dm)

        # --- 2) Run time course with SBML-ID columns ----
        tc = run_time_course(
            start_time=0.0,
            duration=interval,
            intervals=self.intervals,
            update_model=True,
            use_sbml_id=True,   # <-- critical
            model=self.dm,
        )

        # Extract time points
        time = tc["Time"].tolist() if "Time" in tc.columns else []

        # --- 3) Read back final state: export SBML IDs ----
        species_concentrations = {
            sbml_id: _get_transient_concentration(
                name=self.sbml_to_name[sbml_id],
                dm=self.dm
            )
            for sbml_id in self.species_ids
        }

        # --- 4) Reaction fluxes (COPASI reaction IDs already match SBML IDs) ----
        rxn_df = get_reactions(model=self.dm)
        reaction_fluxes = {
            rxn_id: float(rxn_df.loc[rxn_id, "flux"])
            for rxn_id in self.reaction_ids
        }

        return {
            "species_concentrations": species_concentrations,
            "fluxes": reaction_fluxes,
            "time": time,
        }




def run_copasi_utc(core):

    copasi_process = CopasiUTCStep({
                'model_source': 'models/BIOMD0000000012_url.xml',  # represillator model
                'time': 10.0,
                'n_points': 5,
            }, core=core)

    initial_state = copasi_process.initial_state()

    print(f'Initial state: {initial_state}')

    results = copasi_process.update(initial_state)

    print(f'Results: {results}')


def run_copasi_ss(core):

    copasi_process = CopasiSteadyStateStep({
                'model_source': 'models/BIOMD0000000012_url.xml',  # represillator model
            }, core=core)

    initial_state = copasi_process.initial_state()

    print(f'Initial state: {initial_state}')

    results = copasi_process.update(initial_state)

    print(f'Results: {results}')


if __name__ == '__main__':
    from process_bigraph import allocate_core
    core = allocate_core()

    run_copasi_utc(core=core)
    run_copasi_ss(core=core)
