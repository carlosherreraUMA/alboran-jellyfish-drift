"""Gridded ocean-current forcing: regular lon/lat grids, bilinear/linear interpolation,
land masks, bathymetry and multi-level (depth-interpolated) currents."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr
from numba import njit


@njit(cache=True)
def _interp_levels(u_all, v_all, land, depths, lon, lat, z, lon0, lat0, dlon, dlat, ny, nx):  # pragma: no cover
    """Bilinear (horizontal) x linear (vertical) interpolation of stacked levels, one particle at a time.

    ``u_all``/``v_all``/``land`` are flattened (level, lat, lon) stacks on one regular grid. Where the
    lower bracketing level is land (seabed shallower than that level) the upper level is used.
    """
    n = lon.shape[0]
    n_lev = depths.shape[0]
    n_cells = ny * nx
    u = np.empty(n)
    v = np.empty(n)
    for p in range(n):
        fy = (lat[p] - lat0) / dlat
        fx = (lon[p] - lon0) / dlon
        i0 = min(max(int(np.floor(fy)), 0), ny - 2)
        j0 = min(max(int(np.floor(fx)), 0), nx - 2)
        wy = min(max(fy - i0, 0.0), 1.0)
        wx = min(max(fx - j0, 0.0), 1.0)
        base = i0 * nx + j0
        k = 0
        while k + 1 < n_lev and depths[k + 1] <= z[p]:
            k += 1
        k1 = min(k + 1, n_lev - 1)
        frac = 0.0
        if k1 > k:
            frac = min(max((z[p] - depths[k]) / (depths[k1] - depths[k]), 0.0), 1.0)
            nearest = base + (nx if wy >= 0.5 else 0) + (1 if wx >= 0.5 else 0)
            if land[k1 * n_cells + nearest]:
                frac = 0.0
        w00 = (1.0 - wy) * (1.0 - wx)
        w01 = (1.0 - wy) * wx
        w10 = wy * (1.0 - wx)
        w11 = wy * wx
        b0 = k * n_cells + base
        b1 = k1 * n_cells + base
        ua = w00 * u_all[b0] + w01 * u_all[b0 + 1] + w10 * u_all[b0 + nx] + w11 * u_all[b0 + nx + 1]
        va = w00 * v_all[b0] + w01 * v_all[b0 + 1] + w10 * v_all[b0 + nx] + w11 * v_all[b0 + nx + 1]
        ub = w00 * u_all[b1] + w01 * u_all[b1 + 1] + w10 * u_all[b1 + nx] + w11 * u_all[b1 + nx + 1]
        vb = w00 * v_all[b1] + w01 * v_all[b1 + 1] + w10 * v_all[b1 + nx] + w11 * v_all[b1 + nx + 1]
        u[p] = (1.0 - frac) * ua + frac * ub
        v[p] = (1.0 - frac) * va + frac * vb
    return u, v


@dataclass(frozen=True)
class Grid:
    """Regular, ascending lon/lat grid of cell centres."""

    lon: np.ndarray
    lat: np.ndarray

    def __post_init__(self):
        for name, axis in (("lon", self.lon), ("lat", self.lat)):
            steps = np.diff(axis)
            if len(axis) < 2 or np.any(steps <= 0):
                raise ValueError(f"{name} axis must be strictly ascending")
            if np.ptp(steps) > 1e-4 * abs(steps.mean()):
                raise ValueError(f"{name} axis is not regular")

    @property
    def dlon(self) -> float:
        return float((self.lon[-1] - self.lon[0]) / (len(self.lon) - 1))

    @property
    def dlat(self) -> float:
        return float((self.lat[-1] - self.lat[0]) / (len(self.lat) - 1))

    @property
    def shape(self) -> tuple[int, int]:
        return len(self.lat), len(self.lon)

    def fractional_index(self, lon, lat):
        """Continuous (row, col) index of points; integer values are cell centres."""
        return (np.asarray(lat) - self.lat[0]) / self.dlat, (np.asarray(lon) - self.lon[0]) / self.dlon

    def nearest_index(self, lon, lat, weights=None):
        """Indices (row, col) of the cell containing each point, clipped to the grid."""
        i0, j0, wy, wx = weights if weights is not None else self.weights(lon, lat)
        return i0 + (wy >= 0.5), j0 + (wx >= 0.5)

    def contains(self, lon, lat):
        """True for points inside the outer cell centres (where interpolation is defined)."""
        lon = np.asarray(lon)
        lat = np.asarray(lat)
        return (lon >= self.lon[0]) & (lon <= self.lon[-1]) & (lat >= self.lat[0]) & (lat <= self.lat[-1])

    def weights(self, lon, lat):
        """Bilinear stencil (i0, j0, wy, wx) of points, clipped at the grid edges.

        Computing it once and reusing it for u, v and every level with the same grid is the
        main performance lever of the model.
        """
        fy, fx = self.fractional_index(lon, lat)
        ny, nx = self.shape
        i0 = np.minimum(np.maximum(np.floor(fy), 0), ny - 2).astype(np.intp)
        j0 = np.minimum(np.maximum(np.floor(fx), 0), nx - 2).astype(np.intp)
        wy = np.minimum(np.maximum(fy - i0, 0.0), 1.0)
        wx = np.minimum(np.maximum(fx - j0, 0.0), 1.0)
        return i0, j0, wy, wx

    def interpolate(self, field2d: np.ndarray, weights):
        """Bilinear interpolation of a NaN-free 2-D field with a precomputed stencil."""
        i0, j0, wy, wx = weights
        nx = self.shape[1]
        flat = field2d.ravel()
        base = i0 * nx + j0
        f00, f01 = flat[base], flat[base + 1]
        f10, f11 = flat[base + nx], flat[base + nx + 1]
        return (1 - wy) * ((1 - wx) * f00 + wx * f01) + wy * ((1 - wx) * f10 + wx * f11)

    def bilinear(self, field2d: np.ndarray, lon, lat):
        """Bilinear interpolation of a NaN-free 2-D field at points (clipped at the edges)."""
        return self.interpolate(field2d, self.weights(lon, lat))

    def same_as(self, other: Grid) -> bool:
        return self.shape == other.shape and np.allclose(self.lon, other.lon) and np.allclose(self.lat, other.lat)


class GriddedLayer:
    """Time series of (u, v) currents on one vertical level.

    Slices are served from memory or loaded lazily from a sequence of NetCDF files (for
    example one file per month), so multi-year hindcasts fit in a few hundred MB of RAM.
    Land (or seabed shallower than the level) is wherever the first slice is NaN.
    """

    def __init__(self, grid: Grid, times: np.ndarray, loader, land: np.ndarray, name: str = "layer",
                 depth: float = 0.0):
        self.grid = grid
        self.depth = float(depth)  # metres below the surface (0 for the surface layer)
        self.times = np.asarray(times, dtype="datetime64[s]")
        if np.any(np.diff(self.times).astype(np.int64) <= 0):
            raise ValueError("time axis must be strictly ascending")
        self._loader = loader  # k -> (u2d, v2d) with NaNs
        self.land = land
        self.name = name
        self._filled: OrderedDict[int, tuple[np.ndarray, np.ndarray]] = OrderedDict()
        self._at_time: OrderedDict[np.datetime64, tuple[np.ndarray, np.ndarray]] = OrderedDict()

    # ---- construction -------------------------------------------------------------------
    @classmethod
    def from_arrays(cls, times, lon, lat, u, v, name="layer", depth: float = 0.0) -> GriddedLayer:
        u = np.asarray(u, dtype=np.float32)
        v = np.asarray(v, dtype=np.float32)
        grid = Grid(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
        if u.shape != (len(times), *grid.shape) or v.shape != u.shape:
            raise ValueError(f"u/v must have shape (time, lat, lon) = {(len(times), *grid.shape)}")
        land = np.isnan(u[0]) | np.isnan(v[0])
        return cls(grid, times, lambda k: (u[k], v[k]), land, name=name, depth=depth)

    @classmethod
    def from_files(cls, paths: Sequence[str | Path], depth_m: float | None = None, name="layer",
                   cache_files: int = 2) -> GriddedLayer:
        """Open NetCDF files with variables ``uo``/``vo`` (dims time, [depth], latitude, longitude).

        Files are indexed on construction (time axes only) and loaded whole on demand.
        If a depth dimension exists, the level nearest to ``depth_m`` is used and its actual
        depth is stored in ``layer.depth``; files without a depth dimension are surface fields.
        """
        paths = sorted(Path(p) for p in paths)
        if not paths:
            raise FileNotFoundError("no forcing files given")
        index: list[tuple[int, int]] = []  # global k -> (file number, local index)
        all_times: list[np.datetime64] = []
        grid = None
        level = 0.0
        for n, path in enumerate(paths):
            with xr.open_dataset(path) as ds:
                if grid is None:
                    grid = Grid(ds["longitude"].values.astype(np.float64), ds["latitude"].values.astype(np.float64))
                    if "depth" in ds["uo"].dims:
                        d = ds["depth"]
                        level = float(d[0] if depth_m is None else d.sel(depth=depth_m, method="nearest"))
                t = ds["time"].values.astype("datetime64[s]")
            for local, tt in enumerate(t):
                if all_times and tt <= all_times[-1]:
                    continue  # overlapping/duplicate time stamps between files
                all_times.append(tt)
                index.append((n, local))

        cache: OrderedDict[int, tuple[np.ndarray, np.ndarray]] = OrderedDict()

        def load_file(n: int):
            if n in cache:
                cache.move_to_end(n)
                return cache[n]
            with xr.open_dataset(paths[n]) as ds:
                u, v = ds["uo"], ds["vo"]
                if "depth" in u.dims:
                    if depth_m is None:
                        u, v = u.isel(depth=0), v.isel(depth=0)
                    else:
                        u, v = u.sel(depth=depth_m, method="nearest"), v.sel(depth=depth_m, method="nearest")
                arrays = (u.transpose("time", "latitude", "longitude").values.astype(np.float32),
                          v.transpose("time", "latitude", "longitude").values.astype(np.float32))
            cache[n] = arrays
            while len(cache) > cache_files:
                cache.popitem(last=False)
            return arrays

        def loader(k: int):
            n, local = index[k]
            u, v = load_file(n)
            return u[local], v[local]

        u0, v0 = loader(0)
        land = np.isnan(u0) | np.isnan(v0)
        return cls(grid, np.array(all_times), loader, land, name=name, depth=level)

    def cropped_to(self, grid: Grid) -> GriddedLayer | None:
        """This layer restricted to ``grid`` if ``grid`` is an aligned sub-grid of it (same origin and
        spacing, fewer or equal rows/columns); ``None`` otherwise."""
        g = self.grid
        ny, nx = grid.shape
        aligned = (ny <= g.shape[0] and nx <= g.shape[1]
                   and np.allclose(g.lon[:nx], grid.lon, atol=1e-6) and np.allclose(g.lat[:ny], grid.lat, atol=1e-6))
        if not aligned:
            return None
        if g.shape == grid.shape:
            return self
        loader = self._loader

        def cropped(k):
            u, v = loader(k)
            return u[:ny, :nx], v[:ny, :nx]

        return GriddedLayer(grid, self.times, cropped, self.land[:ny, :nx], name=self.name, depth=self.depth)

    # ---- sampling -----------------------------------------------------------------------
    def _slice(self, k: int):
        if k in self._filled:
            self._filled.move_to_end(k)
            return self._filled[k]
        u, v = self._loader(k)
        pair = (np.nan_to_num(u, nan=0.0), np.nan_to_num(v, nan=0.0))
        self._filled[k] = pair
        while len(self._filled) > 4:
            self._filled.popitem(last=False)
        return pair

    def covers(self, t) -> bool:
        t = np.datetime64(t, "s")
        return bool(self.times[0] <= t <= self.times[-1])

    def fields_at(self, t):
        """Full (u, v) grids linearly interpolated to time ``t`` (cached for the RK stages)."""
        t = np.datetime64(t, "s")
        if t in self._at_time:
            return self._at_time[t]
        if not self.covers(t):
            raise ValueError(f"{self.name}: time {t} outside forcing range {self.times[0]}..{self.times[-1]}")
        k = int(np.searchsorted(self.times, t, side="right") - 1)
        k = min(k, len(self.times) - 2)
        t0, t1 = self.times[k], self.times[k + 1]
        w = float((t - t0).astype(np.int64)) / float((t1 - t0).astype(np.int64))
        ua, va = self._slice(k)
        ub, vb = self._slice(k + 1)
        pair = ((1 - w) * ua + w * ub, (1 - w) * va + w * vb)
        self._at_time[t] = pair
        while len(self._at_time) > 6:
            self._at_time.popitem(last=False)
        return pair

    def velocity(self, t, lon, lat, weights=None):
        """(u, v) in m/s at time ``t`` (linear in time, bilinear in space; land counts as 0)."""
        u2d, v2d = self.fields_at(t)
        w = weights if weights is not None else self.grid.weights(lon, lat)
        return self.grid.interpolate(u2d, w), self.grid.interpolate(v2d, w)

    def is_land(self, lon, lat, weights=None):
        """True where the containing cell is land (or seabed above this level)."""
        i, j = self.grid.nearest_index(lon, lat, weights)
        return self.land[i, j]


class Bathymetry:
    """Sea-floor depth (m, positive down) on its own regular grid; land is NaN."""

    def __init__(self, grid: Grid, depth: np.ndarray):
        self.grid = grid
        self.depth = np.asarray(depth, dtype=np.float64)

    @classmethod
    def from_file(cls, path: str | Path, var: str = "deptho") -> Bathymetry:
        with xr.open_dataset(path) as ds:
            d = ds[var].squeeze()
            grid = Grid(ds["longitude"].values.astype(np.float64), ds["latitude"].values.astype(np.float64))
            return cls(grid, d.transpose("latitude", "longitude").values)

    def at(self, lon, lat):
        i, j = self.grid.nearest_index(lon, lat)
        return self.depth[i, j]


class LayeredForcing:
    """Currents at any depth, linearly interpolated between vertical levels.

    ``layers[0]`` must be the surface layer (hourly in practice); deeper layers may have a
    different time step (e.g. daily means). Between two levels the velocity is interpolated
    linearly in depth; below the deepest available level at a location (the level is "land"
    there because the seabed is shallower) the deepest valid level above it is used.
    Stranding is decided by the surface land mask.
    """

    def __init__(self, layers: list[GriddedLayer], bathymetry: Bathymetry | None = None):
        if not layers:
            raise ValueError("at least a surface layer is required")
        layers = sorted(layers, key=lambda layer: layer.depth)
        cropped = [layers[0]] + [layer.cropped_to(layers[0].grid) for layer in layers[1:]]
        self.stacked = all(c is not None for c in cropped)  # all levels on the surface grid
        self.layers = cropped if self.stacked else layers
        self.depths = np.array([layer.depth for layer in self.layers])
        self.bathymetry = bathymetry
        if self.stacked:
            self._land_stack = np.stack([layer.land.ravel() for layer in self.layers])
        self._stack_cache: OrderedDict[np.datetime64, tuple[np.ndarray, np.ndarray]] = OrderedDict()
        # layers sharing a grid share interpolation stencils
        self._grid_of: list[int] = []
        self._grids: list[Grid] = []
        for layer in self.layers:
            match = next((g for g, grid in enumerate(self._grids) if grid.same_as(layer.grid)), None)
            if match is None:
                self._grids.append(layer.grid)
                match = len(self._grids) - 1
            self._grid_of.append(match)

    @property
    def surface(self) -> GriddedLayer:
        return self.layers[0]

    @property
    def grid(self) -> Grid:
        return self.surface.grid

    def velocity(self, t, lon, lat, z=None):
        lon = np.asarray(lon, dtype=np.float64)
        lat = np.asarray(lat, dtype=np.float64)
        if self.stacked:
            z = np.zeros(lon.shape) if z is None else np.broadcast_to(np.asarray(z, dtype=np.float64), lon.shape)
            return self._velocity_stacked(t, lon, lat, z)
        if z is None or len(self.layers) == 1:
            return self.surface.velocity(t, lon, lat)
        z = np.broadcast_to(np.asarray(z, dtype=np.float64), lon.shape)
        stencils = [g.weights(lon, lat) for g in self._grids]

        def sub(n, idx):
            return tuple(a[idx] for a in stencils[self._grid_of[n]])

        u = np.zeros(lon.shape)
        v = np.zeros(lon.shape)
        k = np.searchsorted(self.depths, z, side="right") - 1
        k = np.minimum(np.maximum(k, 0), len(self.depths) - 1)
        for level in range(len(self.layers)):
            idx = np.flatnonzero(k == level)
            if len(idx) == 0:
                continue
            upper = self.layers[level]
            w_up = sub(level, idx)
            ua, va = upper.velocity(t, None, None, w_up)
            if level == len(self.layers) - 1:
                u[idx], v[idx] = ua, va
                continue
            lower = self.layers[level + 1]
            w_low = sub(level + 1, idx)
            frac = np.minimum(np.maximum((z[idx] - upper.depth) / (lower.depth - upper.depth), 0.0), 1.0)
            frac = np.where(lower.is_land(None, None, w_low), 0.0, frac)  # below the local seabed: upper level
            ub, vb = lower.velocity(t, None, None, w_low)
            u[idx] = (1 - frac) * ua + frac * ub
            v[idx] = (1 - frac) * va + frac * vb
        return u, v

    def _stack_at(self, t):
        t = np.datetime64(t, "s")
        if t not in self._stack_cache:
            fields = [layer.fields_at(t) for layer in self.layers]
            self._stack_cache[t] = (np.stack([f[0].ravel() for f in fields]).ravel(),
                                    np.stack([f[1].ravel() for f in fields]).ravel())
            while len(self._stack_cache) > 6:
                self._stack_cache.popitem(last=False)
        return self._stack_cache[t]

    def _velocity_stacked(self, t, lon, lat, z):
        """All levels share one grid: compiled single pass over particles."""
        g = self.grid
        u_all, v_all = self._stack_at(t)
        return _interp_levels(u_all, v_all, self._land_stack.ravel(), self.depths.astype(np.float64),
                              np.ascontiguousarray(lon, dtype=np.float64), np.ascontiguousarray(lat, dtype=np.float64),
                              np.ascontiguousarray(z, dtype=np.float64), float(g.lon[0]), float(g.lat[0]),
                              g.dlon, g.dlat, g.shape[0], g.shape[1])

    def bottom_depth(self, lon, lat):
        """Sea-floor depth (m).

        Without bathymetry, the depth of the shallowest level that is "land" at the location is
        used (infinite if every level is wet) - a crude proxy meant for tests only.
        """
        if self.bathymetry is not None:
            return np.nan_to_num(self.bathymetry.at(lon, lat), nan=0.0)
        bottom = np.full(np.shape(lon), np.inf)
        for layer in reversed(self.layers):
            bottom = np.where(layer.is_land(lon, lat), layer.depth, bottom)
        return bottom

    def is_land(self, lon, lat):
        return self.surface.is_land(lon, lat)

    def covers(self, t) -> bool:
        return all(layer.covers(t) for layer in self.layers)
