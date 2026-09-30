"""Vertical behaviour of virtual jellyfish.

``DielMigration`` implements the diel vertical migration (DVM) of adult *Pelagia noctiluca*
as documented in docs/vertical_behaviour.md:

* night depth ~5 m, day depth of the order of 100-300 m (Berline et al. 2013; Toma et al. 2026);
* migration synchronised with the sun: animals head down after sunrise and up after sunset;
* finite vertical swimming speed (~100 m/h, Berline et al. 2013; El Rahi et al. 2020), so
  the dawn descent and dusk ascent take hours and cross intermediate current layers;
* bathymetric limitation (Berline et al. 2013): animals stay at least ``bottom_clearance_m``
  above the seabed and do not migrate at all where the water is shallower than
  ``no_migration_below_m``;
* optional inter-individual variability of the day depth.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geo import solar_elevation_deg


@dataclass
class FixedDepth:
    """Passive drifter at a constant depth (0 = surface currents)."""

    depth_m: float = 0.0

    def individual_day_depths(self, n: int, rng: np.random.Generator) -> np.ndarray:
        return np.full(n, self.depth_m)

    def target_depth(self, t, lon, lat, day_depth, bottom) -> np.ndarray:
        return np.full(np.shape(lon), self.depth_m)

    def swim(self, z, target, dt_seconds: float) -> np.ndarray:
        return target


@dataclass
class DielMigration:
    day_depth_m: float = 150.0
    night_depth_m: float = 5.0
    speed_m_per_h: float = 100.0
    bottom_clearance_m: float = 100.0
    no_migration_below_m: float = 100.0
    day_depth_rel_sd: float = 0.0
    day_elevation_deg: float = 0.0  # sun above this elevation = day

    def individual_day_depths(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """Preferred day depth of each new individual (normal around the mean, truncated)."""
        if self.day_depth_rel_sd <= 0:
            return np.full(n, self.day_depth_m)
        d = rng.normal(self.day_depth_m, self.day_depth_rel_sd * self.day_depth_m, n)
        return np.clip(d, self.night_depth_m, 2 * self.day_depth_m)

    def is_day(self, t, lon, lat) -> np.ndarray:
        return solar_elevation_deg(np.datetime64(t, "s"), lon, lat) > self.day_elevation_deg

    def target_depth(self, t, lon, lat, day_depth, bottom) -> np.ndarray:
        """Depth each individual is heading for, given light and the local seabed depth."""
        reachable = np.maximum(bottom - self.bottom_clearance_m, self.night_depth_m)
        day_target = np.minimum(day_depth, reachable)
        day_target = np.where(bottom < self.no_migration_below_m, self.night_depth_m, day_target)
        return np.where(self.is_day(t, lon, lat), day_target, self.night_depth_m)

    def swim(self, z, target, dt_seconds: float) -> np.ndarray:
        """Move towards the target depth at the maximum vertical speed."""
        max_dz = self.speed_m_per_h * dt_seconds / 3600.0
        return z + np.clip(target - z, -max_dz, max_dz)


def behaviour_from_config(day_depth_m, cfg: dict | None = None):
    """``None`` -> passive surface drifter; otherwise a DielMigration with config defaults."""
    if day_depth_m is None:
        return FixedDepth(0.0)
    opts = dict(cfg or {})
    return DielMigration(day_depth_m=float(day_depth_m), **opts)
