"""Command-line entry point: ``jellydrift <command> [options]``."""

from __future__ import annotations

import argparse
import itertools
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .behaviour import behaviour_from_config
from .field import Bathymetry, GriddedLayer, LayeredForcing
from .model import DriftParams, Simulation, Source


def load_config(path: str | Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


BEHAVIOUR_OVERRIDES = {"speed_m_per_h": "v", "day_depth_rel_sd": "sd"}


def resolve_spec(cfg: dict, spec: dict) -> dict:
    """Fill a run specification with defaults from the config."""
    out = {"day_depth_m": spec.get("day_depth_m"),
           "horizontal_diffusivity": float(spec.get("horizontal_diffusivity", cfg["drift"]["horizontal_diffusivity"]))}
    out.update({k: spec[k] for k in BEHAVIOUR_OVERRIDES if k in spec})
    return out


def run_name(spec: dict) -> str:
    dd = spec["day_depth_m"]
    name = ("surface" if dd is None else f"dvm{int(dd)}m") + f"_kh{spec['horizontal_diffusivity']:g}"
    for key, tag in BEHAVIOUR_OVERRIDES.items():
        if key in spec:
            name += f"_{tag}{spec[key]:g}"
    return name


def surface_dir(cfg) -> Path:
    return Path(cfg["forcing"]["data_dir"]) / "surface"


def deep_dir(cfg, depth) -> Path:
    return Path(cfg["forcing"]["data_dir"]) / f"deep_{int(depth)}m"


def month_files(directory: Path, start: pd.Timestamp, end: pd.Timestamp) -> list[Path]:
    months = pd.period_range((start - pd.offsets.MonthBegin(1)).to_period("M"), end.to_period("M"), freq="M")
    files = [directory / f"{m.strftime('%Y-%m')}.nc" for m in months]
    return [f for f in files if f.exists()]


# ---- commands ---------------------------------------------------------------------------
def cmd_observations(cfg, args):
    from .observations import fetch_gbif

    o = cfg["observations"]
    out = Path("data/observations") / (o["scientific_name"].lower().replace(" ", "_") + "_gbif.csv")
    df = fetch_gbif(o["scientific_name"], tuple(o["bbox"]), out)
    print(f"{len(df)} records -> {out}")


def cmd_download(cfg, args):
    from .download import download_months

    start = args.start or str((pd.Timestamp(cfg["period"]["start"])
                               - pd.Timedelta(days=cfg["period"]["spinup_days"] + 31)).date())
    end = args.end or cfg["period"]["end"]
    bbox = cfg["domain"]["bbox"]
    if args.what == "surface":
        download_months(cfg["forcing"]["surface_dataset"], surface_dir(cfg), start, end, bbox)
    else:
        depths = [args.depth] if args.depth else cfg["forcing"]["levels_m"]
        for d in depths:
            download_months(cfg["forcing"]["deep_dataset"], deep_dir(cfg, d), start, end, bbox, depth_m=d)


def cmd_check_data(cfg, args):
    """Check that every month is present, readable and that time steps are continuous."""
    from .download import month_starts

    start = str((pd.Timestamp(cfg["period"]["start"]) - pd.Timedelta(days=cfg["period"]["spinup_days"] + 31)).date())
    months = [f"{m:%Y-%m}" for m in month_starts(start, cfg["period"]["end"])]
    dirs = [(surface_dir(cfg), "1h")] + [(deep_dir(cfg, d), "1D") for d in cfg["forcing"]["levels_m"]]
    ok = True
    for directory, step in dirs:
        missing = [m for m in months if not (directory / f"{m}.nc").exists()]
        if missing:
            ok = False
            print(f"{directory}: MISSING {len(missing)} months, e.g. {missing[:3]}")
            continue
        layer = GriddedLayer.from_files([directory / f"{m}.nc" for m in months])
        gaps = np.flatnonzero(np.diff(layer.times) != pd.Timedelta(step).to_timedelta64())
        status = "ok" if len(gaps) == 0 else f"{len(gaps)} GAPS, first after {layer.times[gaps[0]]}"
        ok &= len(gaps) == 0
        print(f"{directory}: {len(months)} months, {layer.times[0]} .. {layer.times[-1]}, "
              f"level {layer.depth:.1f} m, {status}")
    bathy = Path(cfg["forcing"]["bathymetry"])
    print(f"{bathy}: {'ok' if bathy.exists() else 'MISSING'}")
    raise SystemExit(0 if ok and bathy.exists() else 1)


def build_forcing(cfg: dict, start: pd.Timestamp, end: pd.Timestamp, with_depth: bool = True) -> LayeredForcing:
    """Hourly surface currents plus (optionally) the daily 3-D levels and the bathymetry."""
    layers = [GriddedLayer.from_files(month_files(surface_dir(cfg), start, end), name="surface")]
    if with_depth:
        for d in cfg["forcing"]["levels_m"]:
            layers.append(GriddedLayer.from_files(month_files(deep_dir(cfg, d), start, end), depth_m=d,
                                                  name=f"level{d}"))
    bathy = Path(cfg["forcing"]["bathymetry"])
    return LayeredForcing(layers, Bathymetry.from_file(bathy) if bathy.exists() else None)


def run_year(cfg: dict, spec: dict, year: int, experiment: str = "main", out_root: str = "outputs/runs") -> str:
    """Simulate one calendar year (plus spin-up) for one run specification."""
    spec = resolve_spec(cfg, spec)
    day_depth, kh = spec["day_depth_m"], spec["horizontal_diffusivity"]
    out_dir = Path(out_root) / run_name(spec)
    out_file = out_dir / f"strandings_{year}.parquet"
    if out_file.exists():
        return f"{out_file} (exists)"
    p0, p1 = pd.Timestamp(cfg["period"]["start"]), pd.Timestamp(cfg["period"]["end"])
    y0, y1 = max(pd.Timestamp(f"{year}-01-01"), p0), min(pd.Timestamp(f"{year + 1}-01-01"), p1)
    t_start = y0 - pd.Timedelta(days=cfg["period"]["spinup_days"])

    forcing = build_forcing(cfg, t_start, y1, with_depth=day_depth is not None)
    behaviour = behaviour_from_config(day_depth, {**cfg.get("behaviour", {}),
                                                  **{k: spec[k] for k in BEHAVIOUR_OVERRIDES if k in spec}})
    sources = [Source(s["name"], *s["bbox"], particles_per_day=s["particles_per_day"]) for s in cfg["sources"]]
    d = cfg["drift"]
    params = DriftParams(dt_seconds=d["dt_seconds"], horizontal_diffusivity=kh, max_age_days=d["max_age_days"],
                         max_stranding_depth_m=d.get("max_stranding_depth_m", 25.0),
                         seed=int(d["seed"]) * 10_000 + year)
    tic = time.time()
    res = Simulation(forcing, sources, params, behaviour).run(t_start.to_datetime64(), y1.to_datetime64())
    st = res.strandings
    st = st[st["strand_time"] >= y0.to_datetime64()].reset_index(drop=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "spec.json").write_text(json.dumps({"experiment": experiment, **spec}, indent=2))
    st.to_parquet(out_file, index=False)
    res.summary.assign(year=year, run=run_name(spec), seconds=round(time.time() - tic)).to_csv(
        out_dir / f"summary_{year}.csv", index=False)
    return f"{out_file}: {len(st)} strandings in {time.time() - tic:.0f}s"


def cmd_run(cfg, args):
    depth = None if args.day_depth in (None, "none", "surface") else float(args.day_depth)
    spec = {"day_depth_m": depth}
    if args.kh is not None:
        spec["horizontal_diffusivity"] = args.kh
    print(run_year(cfg, spec, args.year, experiment="manual"))


def cmd_sweep(cfg, args):
    p0, p1 = pd.Timestamp(cfg["period"]["start"]), pd.Timestamp(cfg["period"]["end"])
    years = args.years or list(range(p0.year, p1.year + 1))
    names = list(cfg["experiments"]) if args.experiment == "all" else [args.experiment]
    specs = [(name, spec) for name in names for spec in cfg["experiments"][name]]
    # slowest jobs (deepest migration) first so the pool stays busy at the end
    specs.sort(key=lambda x: -(x[1].get("day_depth_m") or 0))
    jobs = [(name, spec, y) for (name, spec), y in itertools.product(specs, years)]
    print(f"{len(jobs)} jobs on {args.workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_year, cfg, spec, y, name): (run_name(resolve_spec(cfg, spec)), y)
                   for name, spec, y in jobs}
        for fut in as_completed(futures):
            try:
                print(f"[done] {futures[fut]} {fut.result()}", flush=True)
            except Exception as exc:  # keep the sweep going, report at the end
                print(f"[FAILED] {futures[fut]} {exc!r}", flush=True)


def cmd_evaluate(cfg, args):
    from .evaluate import evaluate_all

    evaluate_all(cfg, Path(args.runs), Path(args.out))


def cmd_figures(cfg, args):
    from .figures import make_all

    make_all(cfg, Path(args.runs), Path(args.eval), Path(args.out))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="jellydrift", description=__doc__)
    ap.add_argument("--config", default="config/default.yaml")
    sub = ap.add_subparsers(dest="command", required=True)

    sub.add_parser("observations", help="fetch GBIF occurrences")
    sub.add_parser("check-data", help="check downloaded forcing files for missing months and time gaps")
    p = sub.add_parser("download", help="download Copernicus Marine currents (needs login)")
    p.add_argument("what", choices=["surface", "deep"])
    p.add_argument("--depth", type=float)
    p.add_argument("--start")
    p.add_argument("--end")
    p = sub.add_parser("run", help="simulate one year for one parameter set")
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--day-depth", default=None)
    p.add_argument("--kh", type=float, default=None)
    p = sub.add_parser("sweep", help="simulate all runs of an experiment (config 'experiments') in parallel")
    p.add_argument("--experiment", default="main", help="experiment name in the config, or 'all'")
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--years", type=int, nargs="*")
    p = sub.add_parser("evaluate", help="score simulations against observations")
    p.add_argument("--runs", default="outputs/runs")
    p.add_argument("--out", default="outputs/evaluation")
    p = sub.add_parser("figures", help="make report figures")
    p.add_argument("--runs", default="outputs/runs")
    p.add_argument("--eval", default="outputs/evaluation")
    p.add_argument("--out", default="report/figures")

    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    commands = {"observations": cmd_observations, "check-data": cmd_check_data, "download": cmd_download,
                "run": cmd_run, "sweep": cmd_sweep, "evaluate": cmd_evaluate, "figures": cmd_figures}
    commands[args.command](cfg, args)


if __name__ == "__main__":
    main()
