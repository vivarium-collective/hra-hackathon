"""SBML metadata extraction and time unit detection."""

from __future__ import annotations

import libsbml

# Time units known from publications but missing from SBML metadata
_TIME_UNIT_OVERRIDES = {
    "Nyman2011":  "min",
    "Koenig2012": "min",
}


def sbml_time_unit(model):
    """Determine the human-readable time unit from an SBML model."""
    tu_id = model.getTimeUnits() or ""
    if tu_id in ("second", "s"):
        return "s"
    if tu_id in ("minute", "min"):
        return "min"
    if tu_id in ("hour", "h"):
        return "h"
    if tu_id in ("day", "d"):
        return "day"

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

    model_name = (model.getName() or model.getId() or "")
    for key, unit in _TIME_UNIT_OVERRIDES.items():
        if key.lower() in model_name.lower():
            return unit

    return "s"


def get_sbml_info(sbml_path):
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
        "time_unit": sbml_time_unit(m),
    }
