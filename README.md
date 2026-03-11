# HRA Hackathon

Simulation of [Human Reference Atlas](https://humanatlas.io/)-relevant biological models using [COPASI](https://copasi.org/) via [process-bigraph](https://github.com/vivarium-collective/process-bigraph) composites. Models are fetched from [BioModels](https://www.ebi.ac.uk/biomodels/), simulated as time-course experiments, and rendered into an interactive HTML report with species name resolution from UniProt and KEGG.

## Results

**[Download simulation results](https://raw.githubusercontent.com/vivarium-collective/hra-hackathon/main/out_hra_models/results.html)** — interactive Plotly plots with model metadata, hover tooltips, and toggle-able species traces. Save the HTML file and open it in your browser.

| Model ID | Description |
|----------|-------------|
| BIOMD0000000356 | Nyman2011 — Hierarchical insulin-glucose dynamics |
| BIOMD0000000620 | Palmer2014 — IL-1β-blocking therapies in T2DM (disease) |
| BIOMD0000000621 | Palmer2014 — IL-1β-blocking therapies in T2DM (healthy) |
| BIOMD0000000833 | DiCamillo2016 — Insulin signalling pathway |
| BIOMD0000000854 | Gray2016 — The Akt switch model |
| MODEL1209260000 | Koenig2012 — Hepatic glucose metabolism |
| MODEL1912090001 | Huang2014 — Insulin signaling via IRS1 and IRS2 |
| MODEL2401110001 | Pancreas glucose model |

## Setup

Requires Python 3.11–3.12 and [uv](https://docs.astral.sh/uv/).

```sh
uv sync
```

## Usage

Run the default set of HRA models:

```sh
uv run python -m experiments.run_hra_models
```

Run specific models by passing their BioModels IDs:

```sh
uv run python -m experiments.run_hra_models BIOMD0000000356 BIOMD0000000833
```

Or call programmatically:

```python
from process_bigraph import allocate_core
from experiments.run_hra_models import run_hra_models

core = allocate_core()
run_hra_models(core, model_ids=["BIOMD0000000356"])
```

Results are saved to `out_hra_models/results.html`.

## Project structure

```
experiments/
  run_hra_models.py    # Main orchestration — builds documents, runs simulations
  model_ids.py         # Default list of BioModels IDs
  biomodels_fetch.py   # Fetching SBML/SED-ML from BioModels, parsing SED-ML
  species_names.py     # Species name resolution via UniProt & KEGG REST APIs
  sbml_utils.py        # SBML metadata extraction (time units, compartments, etc.)
  plotting.py          # Plotly figure generation and HTML report assembly
processes/
  copasi_process.py    # CopasiUTCStep — process-bigraph Step wrapping COPASI time-course
models/                # Cached SBML files (populated on first run)
out_hra_models/        # Simulation output — JSON state files and results.html
```

## How it works

1. For each model ID, the pipeline fetches metadata and files from BioModels.
2. SED-ML is parsed to extract simulation duration and number of time points.
3. A process-bigraph composite document is built with a `CopasiUTCStep` and a RAM emitter.
4. COPASI runs the time-course simulation (with automatic retry using stiffer solvers if needed).
5. Species names are resolved from SBML annotations via UniProt and KEGG REST APIs.
6. All results are combined into a single interactive HTML report.
