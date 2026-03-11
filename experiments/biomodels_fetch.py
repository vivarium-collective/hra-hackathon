"""BioModels file helpers and SED-ML parsing."""

from __future__ import annotations

import builtins
import os
import re
import tempfile
from pathlib import Path

import biomodels
import libsedml

_SBML_RE = re.compile(r"\.(xml|sbml)$", re.IGNORECASE)
_SEDML_RE = re.compile(r"\.sedml$", re.IGNORECASE)


def file_name(obj):
    return getattr(obj, "name", str(obj))


def iter_entry_files(entry):
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


def fetch_file(file_entry, out_dir):
    f = biomodels.get_file(file_entry)
    if isinstance(f, (str, os.PathLike)) and os.path.exists(str(f)):
        return str(f)
    out_path = os.path.join(out_dir, file_name(file_entry))
    if isinstance(f, bytes):
        Path(out_path).write_bytes(f)
    else:
        Path(out_path).write_text(str(f), encoding="utf-8")
    return out_path


def get_metadata_utf8(model_id):
    """biomodels.get_metadata wrapper -- works around ASCII encoding bug in cache."""
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


def extract_utc(sedml_path):
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


def resolve_sbml_from_sedml(sedml_path, fallback_sbml):
    """If SED-ML references a local SBML file, resolve it; otherwise use fallback."""
    doc = libsedml.readSedMLFromFile(str(sedml_path))
    if doc and doc.getNumModels() > 0:
        src = doc.getModel(0).getSource() or ""
        if src and not src.startswith(("http", "urn:", "biomodels:", "BIOMD")):
            candidate = os.path.join(os.path.dirname(os.path.abspath(sedml_path)), src)
            if os.path.exists(candidate):
                return candidate
    return fallback_sbml


def load_biomodel(model_id, metadata):
    """Fetch SED-ML + SBML from BioModels, return (sbml_path, duration, n_points)."""
    files = list(iter_entry_files(metadata))
    sedml_files = [f for f in files if _SEDML_RE.search(file_name(f))]
    sbml_files = [f for f in files
                  if _SBML_RE.search(file_name(f)) and not _SEDML_RE.search(file_name(f))]

    if not sedml_files or not sbml_files:
        raise ValueError(f"{model_id}: missing SED-ML or SBML in BioModels entry")

    with tempfile.TemporaryDirectory(prefix=f"biomodel_{model_id}_") as tmp:
        sedml_path = fetch_file(sedml_files[0], tmp)
        sbml_path = fetch_file(sbml_files[0], tmp)
        duration, n_points = extract_utc(sedml_path)
        sbml_path = resolve_sbml_from_sedml(sedml_path, sbml_path)

        # Copy to stable location
        stable_dir = os.path.join("models", model_id)
        os.makedirs(stable_dir, exist_ok=True)
        stable_sbml = os.path.join(stable_dir, os.path.basename(sbml_path))
        Path(stable_sbml).write_bytes(Path(sbml_path).read_bytes())

    return stable_sbml, duration, n_points
