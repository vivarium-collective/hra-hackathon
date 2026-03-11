"""Plotting and HTML generation for simulation results."""

from __future__ import annotations

import os
from pathlib import Path

import plotly.graph_objects as go


def extract_result(emitter_data, step_key):
    """Pull the time-course result dict from emitter snapshots."""
    for snapshot in emitter_data:
        if snapshot and step_key in snapshot.get("results", {}):
            r = snapshot["results"][step_key]
            if r:
                return r
    return None


def downsample(time_points, values, max_points=100):
    """Downsample time series to at most max_points, keeping first and last."""
    n = len(time_points)
    if n <= max_points:
        return time_points, values
    step = max(1, n // max_points)
    indices = list(range(0, n, step))
    if indices[-1] != n - 1:
        indices.append(n - 1)
    return [time_points[i] for i in indices], [values[i] for i in indices]


def make_plotly_figure(result_data, title, time_unit="s", name_map=None):
    """Build a Plotly figure from time-course result data."""
    time_pts, vals = downsample(result_data["time"], result_data["values"])
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
