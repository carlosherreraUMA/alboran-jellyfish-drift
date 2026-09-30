"""Public presence-only observations of jellyfish from GBIF (iNaturalist, citizen science, ...)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import requests

from .coast import CoastSegments

GBIF_API = "https://api.gbif.org/v1"

KEEP_COLUMNS = [
    "gbifID", "eventDate", "year", "month", "day", "decimalLongitude", "decimalLatitude",
    "coordinateUncertaintyInMeters", "basisOfRecord", "datasetKey", "datasetName", "license",
    "individualCount", "countryCode",
]


def gbif_predicate(taxon_key: int, bbox: tuple[float, float, float, float]) -> dict:
    """Predicate for an official GBIF download (gives a citable DOI for the exact record set)."""
    lon_min, lon_max, lat_min, lat_max = bbox
    return {
        "type": "and",
        "predicates": [
            {"type": "equals", "key": "TAXON_KEY", "value": str(taxon_key)},
            {"type": "equals", "key": "HAS_COORDINATE", "value": "true"},
            {"type": "equals", "key": "HAS_GEOSPATIAL_ISSUE", "value": "false"},
            {"type": "within", "geometry": (
                f"POLYGON(({lon_min} {lat_min},{lon_max} {lat_min},{lon_max} {lat_max},"
                f"{lon_min} {lat_max},{lon_min} {lat_min}))")},
        ],
    }


def fetch_gbif(scientific_name: str, bbox: tuple[float, float, float, float], out_csv: Path,
               session: requests.Session | None = None) -> pd.DataFrame:
    """Fetch occurrences via the GBIF search API and store a slim CSV (no observer names)."""
    s = session or requests.Session()
    match = s.get(f"{GBIF_API}/species/match", params={"name": scientific_name}, timeout=60).json()
    taxon_key = match["usageKey"]
    lon_min, lon_max, lat_min, lat_max = bbox
    params = {"taxonKey": taxon_key, "decimalLongitude": f"{lon_min},{lon_max}",
              "decimalLatitude": f"{lat_min},{lat_max}", "hasCoordinate": "true",
              "hasGeospatialIssue": "false", "limit": 300, "offset": 0}
    rows = []
    while True:
        page = s.get(f"{GBIF_API}/occurrence/search", params=params, timeout=60).json()
        rows += page["results"]
        if page["endOfRecords"]:
            break
        params["offset"] += params["limit"]
    df = pd.DataFrame(rows)
    df = df[[c for c in KEEP_COLUMNS if c in df.columns]].copy()
    df.insert(0, "scientificName", match.get("scientificName", scientific_name))
    df["occurrenceUrl"] = "https://www.gbif.org/occurrence/" + df["gbifID"].astype(str)
    df = df.sort_values(["eventDate", "gbifID"]).reset_index(drop=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    meta = {"scientificName": scientific_name, "taxonKey": taxon_key, "bbox": bbox,
            "retrieved_utc": pd.Timestamp.now(tz="UTC").isoformat(), "n_records": len(df),
            "download_predicate": gbif_predicate(taxon_key, bbox)}
    out_csv.with_suffix(".json").write_text(json.dumps(meta, indent=2))
    return df


def coastal_events(obs: pd.DataFrame, segments: CoastSegments, max_distance_km: float = 5.0,
                   max_uncertainty_m: float = 5000.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Attach observations to coastal segments and collapse them to (segment, date) events.

    Returns (records with segment/distance/flags, unique events).
    Only records with day-precision dates and coordinate uncertainty <= ``max_uncertainty_m``
    (or unknown) are used.
    """
    df = obs.copy()
    df = df.dropna(subset=["year", "month", "day"])
    df["date"] = pd.to_datetime(dict(year=df["year"], month=df["month"], day=df["day"]), errors="coerce")
    df = df.dropna(subset=["date"])
    unc = df.get("coordinateUncertaintyInMeters")
    if unc is not None:
        df = df[unc.isna() | (unc <= max_uncertainty_m)]
    seg, dist = segments.nearest_segment(df["decimalLongitude"].to_numpy(), df["decimalLatitude"].to_numpy(),
                                         max_distance_km)
    df["segment_id"] = seg
    df["coast_distance_km"] = dist
    matched = df[df["segment_id"] >= 0]
    events = (matched.groupby(["segment_id", "date"]).agg(n_records=("gbifID", "size"))
              .reset_index().sort_values(["date", "segment_id"]).reset_index(drop=True))
    return df.reset_index(drop=True), events
