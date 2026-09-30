"""Monthly downloads of Copernicus Marine current fields (requires `copernicusmarine login`)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def month_starts(start: str, end: str) -> list[pd.Timestamp]:
    return list(pd.date_range(pd.Timestamp(start).to_period("M").to_timestamp(),
                              pd.Timestamp(end).to_period("M").to_timestamp(), freq="MS"))


def download_months(dataset_id: str, out_dir: Path, start: str, end: str, bbox, depth_m: float | None = None,
                    overwrite: bool = False) -> list[Path]:
    """Download ``uo``/``vo`` month by month into ``out_dir/YYYY-MM.nc`` (skips existing files).

    Each file also contains the first time stamp of the following month so consecutive
    files overlap and time interpolation across month boundaries works.
    """
    import copernicusmarine as cm

    lon_min, lon_max, lat_min, lat_max = bbox
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for m0 in month_starts(start, end):
        path = out_dir / f"{m0:%Y-%m}.nc"
        paths.append(path)
        if path.exists() and not overwrite:
            continue
        m1 = m0 + pd.offsets.MonthBegin(1)
        kwargs = dict(
            dataset_id=dataset_id, variables=["uo", "vo"],
            minimum_longitude=lon_min, maximum_longitude=lon_max,
            minimum_latitude=lat_min, maximum_latitude=lat_max,
            start_datetime=f"{m0:%Y-%m-%dT%H:%M:%S}", end_datetime=f"{m1:%Y-%m-%dT%H:%M:%S}",
            output_directory=str(out_dir), output_filename=f"{m0:%Y-%m}.part.nc",
            disable_progress_bar=True,
        )
        if depth_m is not None:
            kwargs.update(minimum_depth=depth_m, maximum_depth=depth_m, coordinates_selection_method="nearest")
        print(f"[download] {dataset_id} {m0:%Y-%m} depth={depth_m} -> {path}", flush=True)
        cm.subset(**kwargs)
        (out_dir / f"{m0:%Y-%m}.part.nc").rename(path)
    return paths
