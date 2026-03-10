# HRA Hackathon

Simulation of HRA-relevant biological models using [COPASI](https://copasi.org/) via [process-bigraph](https://github.com/vivarium-collective/process-bigraph) composites.

## Results

**[View all simulation results](out_hra_models/results.html)** — interactive plots with model details (open the HTML file locally or after cloning).

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

```sh
uv sync
```

## Run

```sh
uv run python -m experiments.run_hra_models
```

Results are saved to `out_hra_models/results.html`.
