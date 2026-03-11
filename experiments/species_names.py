"""Species name resolution via REST APIs (UniProt, KEGG)."""

from __future__ import annotations

import re
import time as _time

import libsbml
import requests

_IDENT_RE = re.compile(r'identifiers\.org/(uniprot|CHEBI|chebi|kegg\.compound)/([A-Za-z0-9:_.-]+)')

# Global cache so we don't re-fetch the same protein across models
_NAME_CACHE: dict[str, str] = {}


def extract_species_identifiers(sbml_path):
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


def fetch_uniprot_name(accession):
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


def fetch_kegg_name(compound_id):
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


def clean_species_name(name):
    """Clean up a species name: strip x##_ prefixes, replace underscores."""
    cleaned = re.sub(r"^x\d+_", "", name)
    cleaned = re.sub(r"\(.*\)$", "", cleaned)
    cleaned = cleaned.replace("_", " ")
    return cleaned.strip() or name


def resolve_species_names(sbml_path):
    """Build a display name mapping for species in an SBML model.

    Returns {copasi_display_name: human_readable_name}.
    """
    doc = libsbml.readSBML(sbml_path)
    model = doc.getModel() if doc else None
    if not model:
        return {}

    identifiers = extract_species_identifiers(sbml_path)
    uniprot_batch = {}

    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        sp_name = sp.getName() or sp.getId()
        sp_ids = identifiers.get(sp_name, {})

        if "uniprot" in sp_ids:
            acc = sp_ids["uniprot"][0]
            if acc not in _NAME_CACHE:
                uniprot_batch[sp_name] = acc

    # Batch fetch UniProt names (with rate limiting)
    for sp_name, acc in uniprot_batch.items():
        fetch_uniprot_name(acc)
        _time.sleep(0.1)

    # Build the final name map
    seen_names: dict[str, int] = {}
    name_map = {}
    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        sp_name = sp.getName() or sp.getId()
        sp_ids = identifiers.get(sp_name, {})

        resolved = ""
        db_label = ""

        if "uniprot" in sp_ids:
            acc = sp_ids["uniprot"][0]
            resolved = _NAME_CACHE.get(acc, "")
            if resolved:
                db_label = f"UniProt:{acc}"

        if not resolved and "kegg.compound" in sp_ids:
            kid = sp_ids["kegg.compound"][0]
            resolved = fetch_kegg_name(kid)
            if resolved:
                db_label = f"KEGG:{kid}"
            _time.sleep(0.1)

        if not resolved:
            sbml_name = sp.getName() or ""
            sbml_id = sp.getId() or ""
            if sbml_name and sbml_name != sbml_id:
                resolved = clean_species_name(sbml_name)
            else:
                resolved = clean_species_name(sbml_id)

        count = seen_names.get(resolved, 0)
        seen_names[resolved] = count + 1
        if count > 0:
            resolved = f"{resolved} ({sp_name})"

        name_map[sp_name] = {"display": resolved, "db_id": db_label}

    return name_map
