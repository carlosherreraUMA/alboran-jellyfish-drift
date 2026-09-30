"""Lagrangian particle-tracking model for jellyfish transport and stranding.

Particles are released continuously from source areas, advected with 4th-order
Runge-Kutta in geographic coordinates, spread with a random walk (horizontal eddy
diffusivity) and removed when they strand on the coast, leave the domain or exceed
their maximum age. Every stranding is recorded, so beaching statistics can be computed
per coastal cell, per day, per source and per particle age.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .behaviour import FixedDepth
from .geo import meters_to_degrees, velocity_to_degrees_per_second


@dataclass
class DriftParams:
    dt_seconds: float = 900.0
    horizontal_diffusivity: float = 1.0  # Kh in m^2/s; random-walk step sigma = sqrt(2 Kh dt)
    drift_factor: float = 1.0  # multiplies the ocean velocity; 1 = fully passive transport
    max_age_days: float = 60.0
    # Only animals in the upper layer can strand; deeper ones meeting the coast (steep slopes)
    # keep their position for that step instead.
    max_stranding_depth_m: float = 25.0
    seed: int = 0


@dataclass
class Source:
    """Rectangular release area (particles are only placed on sea cells)."""

    name: str
    lon_min: float
    lon_max: float
    lat_min: float
    lat_max: float
    particles_per_day: float


@dataclass
class RunResult:
    strandings: pd.DataFrame
    summary: pd.DataFrame
    snapshots: list[pd.DataFrame] = field(default_factory=list)


STRANDED, EXITED, EXPIRED = 1, 2, 3


class Simulation:
    def __init__(self, forcing, sources: list[Source], params: DriftParams | None = None, behaviour=None):
        self.forcing = forcing
        self.sources = sources
        self.params = params or DriftParams()
        self.behaviour = behaviour or FixedDepth(0.0)
        self.rng = np.random.default_rng(self.params.seed)

    # ---- helpers ------------------------------------------------------------------------
    def _rhs(self, t, lon, lat, z=None):
        u, v = self.forcing.velocity(t, lon, lat, z)
        a = self.params.drift_factor
        return velocity_to_degrees_per_second(a * u, a * v, lat)

    def _rk4(self, t, lon, lat, dt_td, z=None):
        """One RK4 step in (lon, lat); depth ``z`` is held fixed during the step."""
        dt = float(dt_td.astype(np.int64))
        half = t + dt_td // 2
        t1 = t + dt_td
        k1x, k1y = self._rhs(t, lon, lat, z)
        k2x, k2y = self._rhs(half, lon + 0.5 * dt * k1x, lat + 0.5 * dt * k1y, z)
        k3x, k3y = self._rhs(half, lon + 0.5 * dt * k2x, lat + 0.5 * dt * k2y, z)
        k4x, k4y = self._rhs(t1, lon + dt * k3x, lat + dt * k3y, z)
        return (lon + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x),
                lat + dt / 6.0 * (k1y + 2 * k2y + 2 * k3y + k4y))

    def integrate(self, lon, lat, start, end, z=None):
        """Pure advection of given positions from ``start`` to ``end`` at fixed depth ``z``
        (no diffusion, no behaviour, no fates).

        Returns the trajectory as arrays of shape (n_steps + 1, n_particles).
        """
        dt_td = np.timedelta64(int(round(self.params.dt_seconds)), "s")
        t = np.datetime64(start, "s")
        lon = np.atleast_1d(np.asarray(lon, dtype=np.float64))
        lat = np.atleast_1d(np.asarray(lat, dtype=np.float64))
        lons, lats = [lon], [lat]
        while t < np.datetime64(end, "s"):
            lon, lat = self._rk4(t, lon, lat, dt_td, z)
            lons.append(lon)
            lats.append(lat)
            t = t + dt_td
        return np.array(lons), np.array(lats)

    def _sample_release(self, src: Source, n: int, max_tries: int = 20):
        lon = np.empty(0)
        lat = np.empty(0)
        for _ in range(max_tries):
            need = n - len(lon)
            if need <= 0:
                break
            cand_lon = self.rng.uniform(src.lon_min, src.lon_max, 2 * need)
            cand_lat = self.rng.uniform(src.lat_min, src.lat_max, 2 * need)
            ok = self.forcing.grid.contains(cand_lon, cand_lat) & ~self.forcing.is_land(cand_lon, cand_lat)
            lon = np.concatenate([lon, cand_lon[ok]])[:n]
            lat = np.concatenate([lat, cand_lat[ok]])[:n]
        return lon, lat

    # ---- main loop ----------------------------------------------------------------------
    def run(self, start, end, release_until=None, snapshot_every_hours: float | None = None,
            snapshot_max_particles: int = 5000, progress: bool = False) -> RunResult:
        p = self.params
        dt = float(p.dt_seconds)
        dt_td = np.timedelta64(int(round(dt)), "s")
        t = np.datetime64(start, "s")
        end = np.datetime64(end, "s")
        release_until = np.datetime64(release_until, "s") if release_until is not None else end
        max_age_s = p.max_age_days * 86400.0
        sigma = np.sqrt(2.0 * p.horizontal_diffusivity * dt) if p.horizontal_diffusivity > 0 else 0.0

        lon = np.empty(0)
        lat = np.empty(0)
        z = np.empty(0)  # depth (m, positive down)
        zday = np.empty(0)  # individual preferred day depth
        born = np.empty(0, dtype="datetime64[s]")
        src_id = np.empty(0, dtype=np.int16)
        pid = np.empty(0, dtype=np.int64)
        next_pid = 0
        carry = np.zeros(len(self.sources))
        released = np.zeros(len(self.sources), dtype=np.int64)
        fates = np.zeros((len(self.sources), 4), dtype=np.int64)
        records: list[dict] = []
        snapshots: list[pd.DataFrame] = []
        next_snapshot = t
        n_steps = int(np.ceil((end - t).astype(np.int64) / dt))

        for step in range(n_steps):
            # 1) release new particles at the start of the step
            if t < release_until:
                for s, src in enumerate(self.sources):
                    carry[s] += src.particles_per_day * dt / 86400.0
                    n_new = int(carry[s])
                    if n_new == 0:
                        continue
                    carry[s] -= n_new
                    new_lon, new_lat = self._sample_release(src, n_new)
                    k = len(new_lon)
                    new_zday = self.behaviour.individual_day_depths(k, self.rng)
                    new_z = self.behaviour.target_depth(t, new_lon, new_lat, new_zday,
                                                        self.forcing.bottom_depth(new_lon, new_lat))
                    lon = np.concatenate([lon, new_lon])
                    lat = np.concatenate([lat, new_lat])
                    z = np.concatenate([z, new_z])
                    zday = np.concatenate([zday, new_zday])
                    born = np.concatenate([born, np.full(k, t)])
                    src_id = np.concatenate([src_id, np.full(k, s, dtype=np.int16)])
                    pid = np.concatenate([pid, np.arange(next_pid, next_pid + k)])
                    next_pid += k
                    released[s] += k

            if snapshot_every_hours and t >= next_snapshot and len(lon):
                sel = slice(None) if len(lon) <= snapshot_max_particles else self.rng.choice(
                    len(lon), snapshot_max_particles, replace=False)
                snapshots.append(pd.DataFrame({"time": t, "pid": pid[sel], "source": src_id[sel],
                                               "lon": lon[sel], "lat": lat[sel], "z": z[sel]}))
                next_snapshot = t + np.timedelta64(int(snapshot_every_hours * 3600), "s")

            if len(lon) == 0:
                t = t + dt_td
                continue

            # 2) vertical behaviour: swim towards the target depth for this step
            target = self.behaviour.target_depth(t, lon, lat, zday, self.forcing.bottom_depth(lon, lat))
            z = self.behaviour.swim(z, target, dt)

            # 3) horizontal advection at that depth: RK4 in (lon, lat)
            t1 = t + dt_td
            new_lon, new_lat = self._rk4(t, lon, lat, dt_td, z)

            # 4) sub-grid dispersion: isotropic random walk
            if sigma > 0:
                dlon, dlat = meters_to_degrees(self.rng.normal(0.0, sigma, len(lon)),
                                               self.rng.normal(0.0, sigma, len(lon)), new_lat)
                new_lon = new_lon + dlon
                new_lat = new_lat + dlat

            # 5) fates
            age_s = (t1 - born).astype(np.int64).astype(np.float64)
            exited = ~self.forcing.grid.contains(new_lon, new_lat)
            hits_land = ~exited & self.forcing.is_land(new_lon, new_lat)
            blocked = hits_land & (z > p.max_stranding_depth_m)
            if np.any(blocked):
                new_lon[blocked] = lon[blocked]
                new_lat[blocked] = lat[blocked]
            stranded = hits_land & ~blocked
            expired = ~exited & ~stranded & (age_s > max_age_s)

            if np.any(stranded):
                i, j = self.forcing.grid.nearest_index(lon[stranded], lat[stranded])
                records.append(pd.DataFrame({
                    "pid": pid[stranded], "source": src_id[stranded],
                    "release_time": born[stranded], "strand_time": t1,
                    "lon": lon[stranded], "lat": lat[stranded],  # last wet position
                    "cell_i": i, "cell_j": j, "age_days": age_s[stranded] / 86400.0, "z": z[stranded],
                }))
            for code, mask in ((STRANDED, stranded), (EXITED, exited), (EXPIRED, expired)):
                if np.any(mask):
                    np.add.at(fates[:, code], src_id[mask], 1)

            keep = ~(stranded | exited | expired)
            lon, lat, z, zday = new_lon[keep], new_lat[keep], z[keep], zday[keep]
            born, src_id, pid = born[keep], src_id[keep], pid[keep]
            t = t1

            if progress and step % 2000 == 0:
                print(f"  step {step}/{n_steps}  t={t}  active={len(lon)}  stranded={fates[:, STRANDED].sum()}",
                      flush=True)

        active = np.bincount(src_id, minlength=len(self.sources))
        summary = pd.DataFrame({
            "source": [s.name for s in self.sources], "released": released,
            "stranded": fates[:, STRANDED], "exited": fates[:, EXITED], "expired": fates[:, EXPIRED],
            "active_at_end": active,
        })
        cols = ["pid", "source", "release_time", "strand_time", "lon", "lat", "cell_i", "cell_j", "age_days", "z"]
        strandings = pd.concat(records, ignore_index=True) if records else pd.DataFrame(columns=cols)
        strandings["source"] = pd.Categorical.from_codes(strandings["source"].astype(int),
                                                         [s.name for s in self.sources])
        return RunResult(strandings=strandings, summary=summary, snapshots=snapshots)
