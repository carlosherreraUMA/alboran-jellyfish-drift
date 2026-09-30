"""Report figures (static PNG) for the stranding hindcast."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from .evaluate import ALL_SOURCES, daily_index, load_events, load_strandings, segments_from_config  # noqa: E402
from .field import GriddedLayer, LayeredForcing  # noqa: E402
from .model import DriftParams, Simulation, Source  # noqa: E402

# Validated reference palette (light mode): categorical slots 1-3, blue sequential ramp, chrome.
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
SOURCE_COLORS = {"offshore_box": SERIES[0], "alboran_basin": SERIES[1], "atlantic_inflow": SERIES[2],
                 ALL_SOURCES: INK2}
SOURCE_LABELS = {"offshore_box": "Offshore box (prototype)", "alboran_basin": "Alboran basin",
                 "atlantic_inflow": "Atlantic inflow", ALL_SOURCES: "All sources pooled"}
BLUES = LinearSegmentedColormap.from_list(
    "ref_blues", ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
LAND = "#e9e8e3"


def _style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "sans-serif", "font.size": 9, "text.color": INK, "axes.labelcolor": INK2,
        "axes.edgecolor": BASELINE, "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "axes.grid": False, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False,
    })


def _save(fig, out: Path, name: str):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[figure] {out / name}.png")


def _surface_layer(cfg, months: list[str]) -> GriddedLayer:
    d = Path(cfg["forcing"]["data_dir"]) / "surface"
    return GriddedLayer.from_files([d / f"{m}.nc" for m in months if (d / f"{m}.nc").exists()], name="surface")


def _map_axes(fig, pos, bbox):
    import cartopy.crs as ccrs

    ax = fig.add_subplot(*pos, projection=ccrs.PlateCarree())
    ax.set_extent(bbox, crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color=GRID)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"color": INK2, "size": 8}
    return ax


def _draw_land(ax, layer: GriddedLayer):
    import cartopy.crs as ccrs

    g = layer.grid
    ax.contourf(g.lon, g.lat, layer.land.astype(float), levels=[0.5, 1.5], colors=[LAND],
                transform=ccrs.PlateCarree(), zorder=1)
    ax.contour(g.lon, g.lat, layer.land.astype(float), levels=[0.5], colors=[BASELINE], linewidths=0.6,
               transform=ccrs.PlateCarree(), zorder=2)


# ---- figures ----------------------------------------------------------------------------
def fig_domain(cfg, segments, out: Path):
    """Model domain, mean surface currents, release areas, coastal segments and observations."""
    import cartopy.crs as ccrs

    layer = _surface_layer(cfg, [f"2022-{m:02d}" for m in range(1, 13)])
    times = layer.times[:-1:6]
    g = layer.grid
    u = np.zeros(g.shape)
    v = np.zeros(g.shape)
    for t in times:  # 6-hourly sample of one year
        uu, vv = layer._slice(int(np.searchsorted(layer.times, t)))
        u += uu
        v += vv
    u, v = u / len(times), v / len(times)
    speed = np.where(layer.land, np.nan, np.hypot(u, v))

    fig = plt.figure(figsize=(10, 5.2))
    ax = _map_axes(fig, (1, 1, 1), cfg["domain"]["bbox"])
    pc = ax.pcolormesh(g.lon, g.lat, speed, cmap=BLUES, vmin=0, vmax=0.8, shading="nearest",
                       transform=ccrs.PlateCarree(), zorder=0)
    _draw_land(ax, layer)
    s = 3
    um, vm = np.where(layer.land, np.nan, u), np.where(layer.land, np.nan, v)
    ax.quiver(g.lon[::s], g.lat[::s], um[::s, ::s], vm[::s, ::s],
              color=INK2, scale=8, width=0.0016, transform=ccrs.PlateCarree(), zorder=3)
    for src in cfg["sources"]:
        x0, x1, y0, y1 = src["bbox"]
        c = SOURCE_COLORS[src["name"]]
        ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], color=c, lw=2, transform=ccrs.PlateCarree(), zorder=4)
        ax.text(x0 + 0.03, y0 + 0.03, SOURCE_LABELS[src["name"]], color=INK, fontsize=8, fontweight="bold",
                transform=ccrs.PlateCarree(), zorder=5,
                bbox=dict(facecolor=SURFACE, edgecolor=c, boxstyle="round,pad=0.2", lw=1))
    tab = segments.table
    for _, r in tab.iterrows():
        if r.segment_id % 3 == 0:
            ax.text(r.lon_c, r.lat_c + 0.09, str(r.segment_id), color=INK2, fontsize=7, ha="center",
                    transform=ccrs.PlateCarree(), zorder=5)
    ax.scatter(tab.lon_c, tab.lat_c, s=9, marker="s", color=INK2, transform=ccrs.PlateCarree(), zorder=5,
               label="Coastal segment (every 3rd numbered)")
    obs = pd.read_csv(Path("data/observations") / (cfg["observations"]["scientific_name"].lower().replace(" ", "_")
                                                   + "_gbif.csv"))
    ax.scatter(obs.decimalLongitude, obs.decimalLatitude, s=10, facecolor="none", edgecolor=INK, lw=0.7,
               transform=ccrs.PlateCarree(), zorder=6, label=f"GBIF P. noctiluca record (n={len(obs)})")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2, fontsize=8)
    cb = fig.colorbar(pc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("Mean surface current speed 2022 (m/s)", color=INK2)
    cb.outline.set_edgecolor(BASELINE)
    ax.set_title("Model domain, mean surface circulation, release areas and observations")
    _save(fig, out, "fig1_domain")


def fig_observations(events: pd.DataFrame, out: Path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 2.8), gridspec_kw={"width_ratios": [1.4, 1]})
    by_year = events.groupby(events.date.dt.year).size()
    by_month = events.groupby(events.date.dt.month).size().reindex(range(1, 13), fill_value=0)
    for ax, series, title in ((axes[0], by_year, "Events per year"), (axes[1], by_month, "Events per month")):
        ax.bar(series.index, series.values, color=SERIES[0], width=0.8)
        ax.grid(axis="y")
        ax.set_axisbelow(True)
        ax.set_title(title)
        ax.set_ylabel("Coastal events (segment-days)")
    axes[1].set_xticks(range(1, 13), list("JFMAMJJASOND"))
    axes[0].set_xticks(by_year.index)
    fig.tight_layout()
    _save(fig, out, "fig2_observations")


RUN_LABELS = {"surface_kh1": "Passive surface drifters", "dvm50m_kh1": "DVM, day depth 50 m",
              "dvm150m_kh1": "DVM, day depth 150 m", "dvm300m_kh1": "DVM, day depth 300 m"}


def _run_label(run: str) -> str:
    return RUN_LABELS.get(run, run)


def fig_skill(skill: pd.DataFrame, age: float, periods: dict, out: Path):
    """Mean percentile score (95% CI) of every run, per source, for calibration and validation."""
    sk = skill[skill.max_age_days == age]
    sources = [s for s in SOURCE_LABELS if s in set(sk.source)]

    def order(r):
        head = r.split("_")[0]
        return (head != "surface", float("".join(c for c in head if c.isdigit()) or 0), r)

    runs = sorted(sk.run.unique(), key=order)
    fig, axes = plt.subplots(1, len(sources), figsize=(2.6 * len(sources) + 1.2, 3.6), sharey=True, sharex=True)
    for ax, src in zip(np.atleast_1d(axes), sources, strict=True):
        for period, color, dy in (("calibration", SERIES[0], -0.15), ("validation", SERIES[1], 0.15)):
            d = sk[(sk.source == src) & (sk.period == period)].set_index("run").reindex(runs)
            y = np.arange(len(runs)) + dy
            y0, y1 = periods[period]
            ax.errorbar(d.mean_score, y, xerr=[d.mean_score - d.ci_low, d.ci_high - d.mean_score], fmt="o",
                        ms=4.5, color=color, ecolor=color, elinewidth=1.2, capsize=0,
                        label=f"{period.capitalize()} ({y0}-{y1})")
        ax.axvline(0.5, color=MUTED, lw=1, ls="--")
        ax.set_title(SOURCE_LABELS[src], fontsize=9)
        ax.grid(axis="x")
        ax.set_axisbelow(True)
        ax.set_xlabel("Mean percentile score")
    ax0 = np.atleast_1d(axes)[0]
    ax0.set_yticks(range(len(runs)), runs)
    ax0.invert_yaxis()
    handles, labels = ax0.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.08))
    fig.suptitle(f"Skill against GBIF events; 0.5 = no skill (dashed); particles <= {age:g} days; 95% bootstrap CI",
                 y=1.13, fontsize=9, color=INK2)
    fig.tight_layout()
    _save(fig, out, "fig3_skill")


def fig_case_study(segments, indices: dict[str, pd.DataFrame], events: pd.DataFrame, out: Path,
                   start="2022-03-10", end="2022-04-10", name="fig4_case_march2022"):
    """Small multiples of the stranding index (one per behaviour) around the March 2022 press report."""
    press = pd.read_csv("data/observations/press_events.csv", parse_dates=["date"])
    seg_press, _ = segments.nearest_segment(press.lon.to_numpy(), press.lat.to_numpy(), 10.0)
    keep = [s for s in segments.table.segment_id if s <= 16]  # Strait of Gibraltar to Salobrena
    labels = segments.table.set_index("segment_id").loc[keep, "label"]
    windows = {run: idx.loc[start:end, keep] for run, idx in indices.items()}
    vmax = max(max(np.nanpercentile(w.to_numpy(), 99) for w in windows.values()), 1e-6)
    ev = events[(events.date >= start) & (events.date <= end) & events.segment_id.isin(keep)]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.4), sharex=True, sharey=True)
    for ax, (run, w) in zip(axes.ravel(), windows.items(), strict=False):
        im = ax.imshow(w.T.to_numpy(), aspect="auto", cmap=BLUES, vmin=0, vmax=vmax, interpolation="nearest",
                       extent=(-0.5, len(w) - 0.5, len(keep) - 0.5, -0.5))
        ax.scatter([(d - w.index[0]).days for d in ev.date], [keep.index(s) for s in ev.segment_id], marker="o",
                   s=40, facecolor="none", edgecolor=SERIES[1], lw=1.6, label="GBIF observation", zorder=3)
        for d, s in zip(press.date, seg_press, strict=True):
            if s in keep and pd.Timestamp(start) <= d <= pd.Timestamp(end):
                ax.scatter([(d - w.index[0]).days], [keep.index(s)], marker="D", s=46, facecolor="none",
                           edgecolor=INK, lw=1.6, zorder=4, label="Press report (Fuengirola, 29 Mar 2022)")
        ax.set_title(_run_label(run), fontsize=9)
        ticks = np.arange(0, len(w), 7)
        ax.set_xticks(ticks, [d.strftime("%d %b") for d in w.index[ticks]])
        ax.set_yticks(range(len(keep)), labels, fontsize=7)
    for ax in axes.ravel()[len(windows):]:
        ax.set_visible(False)
    cb = fig.colorbar(im, ax=axes, shrink=0.6, pad=0.02)
    cb.set_label("Strandings per 1000 particles released/day (all sources, age <= 40 d)", color=INK2)
    cb.outline.set_edgecolor(BASELINE)
    handles, lab = axes.ravel()[0].get_legend_handles_labels()
    uniq = dict(zip(lab, handles, strict=True))
    fig.legend(uniq.values(), uniq.keys(), loc="lower center", ncol=2, bbox_to_anchor=(0.45, -0.02))
    fig.suptitle("Hindcast stranding index on the western Alboran coast, 10 March - 10 April 2022", x=0.45,
                 fontsize=10, fontweight="bold")
    _save(fig, out, name)


def fig_segment_timeseries(indices: dict[str, pd.DataFrame], events: pd.DataFrame, segments, segment_id: int,
                           year: int, out: Path):
    ev = events[(events.segment_id == segment_id) & (events.date.dt.year == year)]
    fig, ax = plt.subplots(figsize=(9, 3))
    top = 0.0
    for color, (run, idx) in zip(SERIES, indices.items(), strict=False):
        s = idx.loc[f"{year}-01-01":f"{year}-12-31", segment_id]
        ax.plot(s.index, s.values, color=color, lw=1.2, drawstyle="steps-mid", label=_run_label(run))
        top = max(top, float(s.max()))
    top = top if top > 0 else 1.0
    ax.vlines(ev.date, 0, top * 1.05, color=INK, lw=1.4, label="Observed (GBIF)")
    ax.set_ylim(0, top * 1.1)
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    ax.set_ylabel("Strandings per 1000\nparticles released/day")
    label = segments.table.set_index("segment_id").loc[segment_id, "label"]
    ax.set_title(f"Segment {label}, {year}: modelled stranding index (all sources, age <= 40 d) and observed events")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8)
    fig.tight_layout()
    _save(fig, out, f"fig5_timeseries_seg{segment_id:02d}_{year}")


def fig_seasonality(indices: dict[str, pd.DataFrame], events: pd.DataFrame, segments, out: Path):
    """Mean monthly stranding index per segment for two behaviours, next to observed events."""
    panels = []
    for run, idx in indices.items():
        panels.append((_run_label(run), idx.groupby(idx.index.month).mean().T, "Strandings per 1000 particles/day"))
    first = panels[0][1]
    obs = events.groupby([events.segment_id, events.date.dt.month]).size().unstack(fill_value=0)
    obs = obs.reindex(index=first.index, columns=range(1, 13), fill_value=0)
    panels.append(("Observed: GBIF events", obs, "Events (segment-days)"))
    labels = segments.table.set_index("segment_id").loc[first.index, "label"]
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 6.8), sharey=True)
    for ax, (title, data, cbl) in zip(axes, panels, strict=True):
        im = ax.imshow(data.to_numpy(), aspect="auto", cmap=BLUES, vmin=0, interpolation="nearest")
        ax.set_xticks(range(12), list("JFMAMJJASOND"))
        ax.set_title(title, fontsize=9)
        cb = fig.colorbar(im, ax=ax, shrink=0.45, pad=0.07, location="bottom")
        cb.set_label(cbl, color=INK2, fontsize=8)
        cb.outline.set_edgecolor(BASELINE)
    axes[0].set_yticks(range(len(labels)), labels, fontsize=7)
    fig.suptitle("Seasonal cycle by coastal segment (model: all sources, age <= 40 d, 2017-2026)", fontsize=10,
                 fontweight="bold")
    _save(fig, out, "fig6_seasonality")


def fig_snapshot(cfg, out: Path, start="2022-03-01", end="2022-03-29"):
    """Particle positions coloured by source after a short hindcast (illustration)."""
    import cartopy.crs as ccrs

    months = sorted({f"{d:%Y-%m}" for d in pd.date_range(start, end, freq="D")})
    layer = _surface_layer(cfg, months)
    sources = [Source(s["name"], *s["bbox"], particles_per_day=s["particles_per_day"]) for s in cfg["sources"]]
    sim = Simulation(LayeredForcing([layer]), sources, DriftParams(dt_seconds=1800, horizontal_diffusivity=1.0,
                                                              max_age_days=40, seed=7))
    res = sim.run(np.datetime64(start), np.datetime64(end), snapshot_every_hours=24, snapshot_max_particles=20000)
    snap = res.snapshots[-1]
    t = pd.Timestamp(snap.time.iloc[0])
    fig = plt.figure(figsize=(10, 4.8))
    ax = _map_axes(fig, (1, 1, 1), cfg["domain"]["bbox"])
    _draw_land(ax, layer)
    g = layer.grid
    k = int(np.searchsorted(layer.times, np.datetime64(t)))
    uu, vv = layer._slice(k)
    s = 3
    ax.quiver(g.lon[::s], g.lat[::s], np.where(layer.land, np.nan, uu)[::s, ::s],
              np.where(layer.land, np.nan, vv)[::s, ::s], color=MUTED, scale=10, width=0.0014,
              transform=ccrs.PlateCarree(), zorder=2)
    names = [s.name for s in sources]
    for sid, name in enumerate(names):
        p = snap[snap.source == sid]
        ax.scatter(p.lon, p.lat, s=2, color=SOURCE_COLORS[name], alpha=0.7, lw=0, transform=ccrs.PlateCarree(),
                   zorder=3, label=f"{SOURCE_LABELS[name]} (n={len(p)})")
    st = res.strandings
    ax.scatter(st.lon, st.lat, s=3, color=INK, alpha=0.25, lw=0, transform=ccrs.PlateCarree(), zorder=4,
               label=f"Last wet position of stranded particles since {start} (n={len(st)})")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2, markerscale=5, fontsize=8)
    ax.set_title(f"Simulated particles on {t:%d %b %Y %H:%M} UTC (surface drifters, 28-day hindcast)")
    _save(fig, out, "fig7_snapshot_march2022")


def fig_vertical_behaviour(cfg, out: Path, day_depth=150.0, start="2022-03-01", end="2022-03-15"):
    """Simulated depth of individual jellyfish (diel vertical migration) and its daily cycle."""
    from .behaviour import behaviour_from_config
    from .cli import build_forcing
    from .geo import solar_elevation_deg

    t0, t1 = pd.Timestamp(start), pd.Timestamp(end)
    forcing = build_forcing(cfg, t0, t1, with_depth=True)
    sources = [Source(s["name"], *s["bbox"], particles_per_day=s["particles_per_day"]) for s in cfg["sources"]]
    sim = Simulation(forcing, sources, DriftParams(dt_seconds=cfg["drift"]["dt_seconds"], horizontal_diffusivity=1.0,
                                                   max_age_days=40, seed=11),
                     behaviour_from_config(day_depth, cfg["behaviour"]))
    res = sim.run(t0.to_datetime64(), t1.to_datetime64(), snapshot_every_hours=0.5, snapshot_max_particles=10**6)
    snaps = pd.concat(res.snapshots, ignore_index=True)
    snaps["time"] = pd.to_datetime(snaps["time"])
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [1.7, 1]})

    ax = axes[0]
    w0, w1 = t1 - pd.Timedelta(days=3), t1
    win = snaps[(snaps.time >= w0) & (snaps.time < w1)]
    alive = win.groupby("pid").time.agg(["min", "max"])
    full = alive[(alive["min"] <= w0) & (alive["max"] >= w1 - pd.Timedelta(hours=1))].index
    rng = np.random.default_rng(0)
    for pid in rng.choice(full, size=min(12, len(full)), replace=False):
        d = win[win.pid == pid]
        ax.plot(d.time, d.z, color=SERIES[0], lw=1, alpha=0.6)
    hours = pd.date_range(w0, w1, freq="10min")
    night = solar_elevation_deg(hours.values.astype("datetime64[s]"), -4.0, 36.2) <= 0
    ax.fill_between(hours, 0, 1, where=night, color=GRID, alpha=0.8, lw=0, transform=ax.get_xaxis_transform(),
                    label="Night (sun below horizon)")
    ax.set_ylim(260, 0)
    ax.set_ylabel("Depth (m)")
    ax.set_title(f"Twelve simulated individuals, mean day depth {day_depth:g} m")
    ax.legend(loc="lower left", fontsize=8)
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%d %b\n%H:%M"))

    ax = axes[1]
    snaps["hour"] = snaps.time.dt.hour + snaps.time.dt.minute / 60
    q = snaps.groupby("hour").z.quantile([0.1, 0.5, 0.9]).unstack()
    ax.fill_between(q.index, q[0.1], q[0.9], color="#cde2fb", lw=0, label="10th-90th percentile")
    ax.plot(q.index, q[0.5], color=SERIES[0], lw=2, label="Median")
    ax.set_ylim(260, 0)
    ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 6))
    ax.set_xlabel("Hour (UTC)")
    ax.set_title("Depth by time of day (all individuals)")
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    ax.legend(loc="lower center", fontsize=8)
    fig.tight_layout()
    _save(fig, out, "fig8_vertical_behaviour")


def make_all(cfg, runs_dir: Path, eval_dir: Path, out: Path):
    _style()
    ev = cfg["evaluation"]
    fcfg = cfg["figures"]
    segments = segments_from_config(cfg)
    events = load_events(cfg, segments)
    fig_domain(cfg, segments, out)
    fig_observations(events, out)
    fig_snapshot(cfg, out)
    fig_vertical_behaviour(cfg, out)
    skill_file = eval_dir / "skill.csv"
    if not skill_file.exists():
        print("[figures] no evaluation yet: skipping skill-dependent figures")
        return
    periods = {"calibration": tuple(ev["calibration_years"]), "validation": tuple(ev["validation_years"])}
    fig_skill(pd.read_csv(skill_file), fcfg["display_age_days"], periods, out)
    dates = pd.date_range(cfg["period"]["start"], cfg["period"]["end"], freq="D")
    ppd = {s["name"]: s["particles_per_day"] for s in cfg["sources"]}
    indices = {run: daily_index(load_strandings(runs_dir / run), segments, dates, ppd,
                                fcfg["display_age_days"])[fcfg["display_source"]]
               for run in fcfg["display_runs"] if (runs_dir / run).exists()}
    fig_case_study(segments, indices, events, out)
    pair = {r: indices[r] for r in fcfg["timeseries_runs"] if r in indices}
    busiest = events.groupby("segment_id").size().idxmax()
    year = int(events[events.segment_id == busiest].date.dt.year.value_counts().idxmax())
    fig_segment_timeseries(pair, events, segments, int(busiest), year, out)
    fig_seasonality(pair, events, segments, out)
