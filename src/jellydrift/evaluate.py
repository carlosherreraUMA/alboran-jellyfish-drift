"""Skill of simulated stranding indices against presence-only observations.

For each observed event (coastal segment s, date d) the modelled stranding index at (s, d)
is ranked against the index at the same segment on the surrounding days d +- W. With no
skill, the (mid-)rank percentile is ~0.5 whatever the seasonality of observer effort, because
events are only compared with nearby days at the same place. Significance comes from a
permutation test that moves each event to a random day of its own window.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .coast import CoastSegments, build_segments
from .field import GriddedLayer
from .observations import coastal_events

ALL_SOURCES = "all_sources"


def segments_from_config(cfg) -> CoastSegments:
    first = sorted((Path(cfg["forcing"]["data_dir"]) / "surface").glob("*.nc"))[0]
    layer = GriddedLayer.from_files([first])
    c = cfg["coast"]
    return build_segments(layer.grid, layer.land, c["lat_min"], c["lon_min"], c["lon_max"], c["segment_width_deg"])


def load_events(cfg, segments: CoastSegments) -> pd.DataFrame:
    o = cfg["observations"]
    obs = pd.read_csv(Path("data/observations") / (o["scientific_name"].lower().replace(" ", "_") + "_gbif.csv"))
    _, events = coastal_events(obs, segments, o["max_coast_distance_km"], o["max_uncertainty_m"])
    p0, p1 = pd.Timestamp(cfg["period"]["start"]), pd.Timestamp(cfg["period"]["end"])
    return events[(events["date"] >= p0) & (events["date"] <= p1)].reset_index(drop=True)


def load_strandings(run_dir: Path) -> pd.DataFrame:
    files = sorted(run_dir.glob("strandings_*.parquet"))
    return pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)


def daily_index(strandings: pd.DataFrame, segments: CoastSegments, dates: pd.DatetimeIndex,
                particles_per_day: dict[str, float], max_age_days: float) -> dict[str, pd.DataFrame]:
    """Strandings per day and segment per 1000 particles released per day, for each source.

    Returns {source: DataFrame(index=dates, columns=segment ids)}, plus the pooled
    ``all_sources`` index normalised by the total release rate.
    """
    st = strandings[strandings["age_days"] <= max_age_days]
    seg = segments.segment_of_cells(st["cell_i"].to_numpy(), st["cell_j"].to_numpy())
    st = st.assign(segment_id=seg, date=pd.to_datetime(st["strand_time"]).dt.floor("D"))
    st = st[st["segment_id"] >= 0]
    n_seg = len(segments.table)
    out = {}
    counts = st.groupby(["source", "date", "segment_id"], observed=True).size()
    for src, ppd in particles_per_day.items():
        c = counts.xs(src, level="source") if src in counts.index.get_level_values(0) else pd.Series(dtype=float)
        mat = c.unstack("segment_id", fill_value=0) if len(c) else pd.DataFrame(index=dates)
        out[src] = mat.reindex(index=dates, columns=range(n_seg), fill_value=0).astype(float) * 1000.0 / ppd
    total = sum(particles_per_day.values())
    out[ALL_SOURCES] = sum(out[s] * particles_per_day[s] for s in particles_per_day) / total
    return out


def reference_matrix(index: pd.DataFrame, events: pd.DataFrame, window: int, lag: int) -> np.ndarray:
    """Rows: events; columns: days d-W..d+W (centre = event day). NaN outside the simulated period.

    The index is first smoothed with a trailing max over ``lag`` days, so a stranding up to
    ``lag`` days before an observation counts; the same smoothing applies to reference days.
    """
    smooth = index.rolling(lag + 1, min_periods=1).max() if lag > 0 else index
    values = smooth.to_numpy()
    day0 = smooth.index[0]
    n_days = len(smooth.index)
    offsets = np.arange(-window, window + 1)
    rows = ((events["date"] - day0).dt.days.to_numpy()[:, None] + offsets[None, :])
    valid = (rows >= 0) & (rows < n_days)
    ref = np.full(rows.shape, np.nan)
    seg = np.broadcast_to(events["segment_id"].to_numpy()[:, None], rows.shape)
    ref[valid] = values[rows[valid], seg[valid]]
    return ref


def midrank_percentiles(ref: np.ndarray) -> np.ndarray:
    """Percentile of every entry within its row (ties count half), NaN-aware."""
    x = ref[:, :, None]
    y = ref[:, None, :]
    with np.errstate(invalid="ignore"):
        less = np.nansum(y < x, axis=2)
        equal = np.nansum(y == x, axis=2)
    n = np.sum(~np.isnan(ref), axis=1, keepdims=True)
    pct = (less + 0.5 * equal) / n
    return np.where(np.isnan(ref), np.nan, pct)


def score_events(ref: np.ndarray, n_permutations: int = 10000, n_bootstrap: int = 2000, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    centre = ref.shape[1] // 2
    ok = ~np.isnan(ref[:, centre]) & (np.sum(~np.isnan(ref), axis=1) > centre)
    ref = ref[ok]
    n = len(ref)
    if n == 0:
        return dict(n_events=0)
    pct = midrank_percentiles(ref)
    scores = pct[:, centre]
    mean = float(scores.mean())
    # permutation null: each event moved to a random valid day of its own window
    valid_cols = [np.flatnonzero(~np.isnan(r)) for r in ref]
    picks = np.stack([rng.choice(v, n_permutations) for v in valid_cols])  # (n_events, n_perm)
    null = pct[np.arange(n)[:, None], picks].mean(axis=0)
    boot = scores[rng.integers(0, n, (n_bootstrap, n))].mean(axis=1)
    hits = ref[:, centre] > 0
    base = np.nanmean(ref > 0, axis=1)
    return dict(n_events=n, mean_score=mean, ci_low=float(np.quantile(boot, 0.025)),
                ci_high=float(np.quantile(boot, 0.975)),
                p_value=float((1 + np.sum(null >= mean)) / (1 + n_permutations)),
                hit_rate=float(hits.mean()), base_rate=float(base.mean()))


def evaluate_all(cfg: dict, runs_dir: Path, out_dir: Path) -> pd.DataFrame:
    ev_cfg = cfg["evaluation"]
    segments = segments_from_config(cfg)
    events = load_events(cfg, segments)
    dates = pd.date_range(cfg["period"]["start"], cfg["period"]["end"], freq="D")
    ppd = {s["name"]: s["particles_per_day"] for s in cfg["sources"]}
    periods = {
        "calibration": tuple(ev_cfg["calibration_years"]),
        "validation": tuple(ev_cfg["validation_years"]),
        "all": (dates[0].year, dates[-1].year),
    }
    rows = []
    for run_dir in sorted(p for p in runs_dir.iterdir() if p.is_dir()):
        spec_file = run_dir / "spec.json"
        experiment = json.loads(spec_file.read_text()).get("experiment", "main") if spec_file.exists() else "main"
        strandings = load_strandings(run_dir)
        for max_age in ev_cfg["max_age_days"]:
            indices = daily_index(strandings, segments, dates, ppd, max_age)
            for src, idx in indices.items():
                for period, (y0, y1) in periods.items():
                    sel = events[(events["date"].dt.year >= y0) & (events["date"].dt.year <= y1)]
                    ref = reference_matrix(idx, sel, ev_cfg["window_days"], ev_cfg["lag_days"])
                    res = score_events(ref, ev_cfg["n_permutations"])
                    rows.append(dict(run=run_dir.name, experiment=experiment, source=src, max_age_days=max_age,
                                     period=period, **res))
        print(f"[evaluate] {run_dir.name} done", flush=True)
    table = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_dir / "skill.csv", index=False)
    events.merge(segments.table[["segment_id", "label"]]).to_csv(out_dir / "events.csv", index=False)
    segments.table.to_csv(out_dir / "segments.csv", index=False)

    # Model selection among the main hypotheses on the calibration period only; report the
    # out-of-sample validation score of the selected configuration.
    key = ["run", "source", "max_age_days"]
    main = table[table["experiment"] == "main"]
    cal = main[main["period"] == "calibration"].sort_values("mean_score", ascending=False)
    best = cal.iloc[0][key].to_dict()
    val = main[(main["period"] == "validation")].set_index(key).loc[tuple(best.values())].drop(
        "experiment").to_dict()
    summary = {"selected_on_calibration": best,
               "calibration": cal.iloc[0].drop(key + ["period", "experiment"]).to_dict(),
               "validation": val,
               "n_configurations": int(len(cal))}
    (out_dir / "selection.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps(summary, indent=2, default=float))
    return table
