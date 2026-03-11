import os
from pathlib import Path
from typing import Dict, Any

from pandas import DataFrame
from process_bigraph import Step
import COPASI
from basico import (
    load_model,
    get_species,
    get_reactions,
    run_time_course,
)

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

        # Retry with stiffer solvers if LSODA produces NaN
        if tc.isnull().any().any():
            print("  LSODA produced NaN, retrying with more internal steps...")
            tc = run_time_course(**{**tc_kwargs, 'max_steps': 500000})

        if tc.isnull().any().any():
            print("  Still NaN, retrying with RADAU5 (stiff solver)...")
            tc = run_time_course(**{**tc_kwargs, 'method': 'radau5'})

        if tc.isnull().any().any():
            print("  Still NaN, retrying RADAU5 with looser tolerances...")
            tc = run_time_course(**{
                **tc_kwargs, 'method': 'radau5',
                'r_tol': 1e-3, 'a_tol': 1e-6, 'max_steps': 10000000,
            })

        time_list = tc.index.to_list()

        result = {
            "time": time_list,
            "columns": [self.sbml_to_name.get(c, c) for c in tc.columns],
            "values": tc.values.tolist(),
        }

        return {"result": result}
