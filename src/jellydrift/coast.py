"""Coastal cells of the ocean-model grid and their aggregation into named coastal segments."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .field import Grid
from .geo import haversine_m

# Approximate positions used only to give segments human-readable labels.
REFERENCE_TOWNS: list[tuple[str, float, float]] = [
    ("Tarifa", -5.60, 36.01),
    ("Algeciras", -5.45, 36.13),
    ("Gibraltar", -5.35, 36.14),
    ("La Linea", -5.35, 36.17),
    ("Sotogrande", -5.27, 36.29),
    ("Manilva", -5.23, 36.37),
    ("Estepona", -5.15, 36.43),
    ("San Pedro de Alcantara", -4.99, 36.48),
    ("Marbella", -4.88, 36.51),
    ("Calahonda", -4.72, 36.49),
    ("Mijas Costa", -4.68, 36.50),
    ("Fuengirola", -4.62, 36.54),
    ("Benalmadena", -4.56, 36.60),
    ("Torremolinos", -4.50, 36.62),
    ("Malaga", -4.42, 36.72),
    ("Rincon de la Victoria", -4.28, 36.72),
    ("Torre del Mar", -4.09, 36.74),
    ("Torrox Costa", -3.95, 36.73),
    ("Nerja", -3.88, 36.75),
    ("Almunecar", -3.69, 36.73),
    ("Salobrena", -3.59, 36.74),
    ("Motril", -3.52, 36.72),
    ("Castell de Ferro", -3.36, 36.72),
    ("La Rabita", -3.17, 36.75),
    ("Adra", -3.02, 36.75),
    ("Balerma", -2.88, 36.73),
    ("Almerimar", -2.80, 36.70),
    ("Roquetas de Mar", -2.61, 36.76),
    ("Almeria", -2.46, 36.84),
    ("Cabo de Gata", -2.19, 36.76),
    ("San Jose", -2.11, 36.76),
    ("Carboneras", -1.90, 36.99),
]


def coastal_cells(land: np.ndarray) -> np.ndarray:
    """Sea cells that touch land (8-neighbourhood)."""
    padded = np.pad(land, 1, constant_values=False)
    touch = np.zeros_like(land)
    ny, nx = land.shape
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                touch |= padded[1 + di:1 + di + ny, 1 + dj:1 + dj + nx]
    return ~land & touch


@dataclass
class CoastSegments:
    table: pd.DataFrame  # segment_id, label, lon_min, lon_max, lon_c, lat_c, n_cells
    cell_segment: np.ndarray  # (ny, nx) int, -1 outside any segment
    grid: Grid

    def segment_of_cells(self, i, j):
        return self.cell_segment[np.asarray(i), np.asarray(j)]

    def nearest_segment(self, lon, lat, max_distance_km: float):
        """Segment of the nearest coastal cell within ``max_distance_km`` (else -1) and the distance."""
        ii, jj = np.nonzero(self.cell_segment >= 0)
        clon, clat = self.grid.lon[jj], self.grid.lat[ii]
        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)
        d = haversine_m(lon[:, None], lat[:, None], clon[None, :], clat[None, :]) / 1000.0
        k = np.argmin(d, axis=1)
        dist = d[np.arange(len(lon)), k]
        seg = np.where(dist <= max_distance_km, self.cell_segment[ii[k], jj[k]], -1)
        return seg, dist


def build_segments(grid: Grid, land: np.ndarray, lat_min: float = 36.0, lon_min: float = -5.65,
                   lon_max: float = -1.5, width_deg: float = 0.125) -> CoastSegments:
    """Group coastal cells north of ``lat_min`` (the Iberian coast) into longitude bins."""
    coast = coastal_cells(land)
    lat2d = np.broadcast_to(grid.lat[:, None], land.shape)
    lon2d = np.broadcast_to(grid.lon[None, :], land.shape)
    selected = coast & (lat2d >= lat_min) & (lon2d >= lon_min) & (lon2d < lon_max)
    edges = np.arange(lon_min, lon_max + width_deg, width_deg)
    bins = np.digitize(lon2d, edges) - 1
    cell_segment = np.full(land.shape, -1, dtype=np.int64)
    rows = []
    for b in range(len(edges) - 1):
        m = selected & (bins == b)
        if not m.any():
            continue
        sid = len(rows)
        cell_segment[m] = sid
        lon_c, lat_c = float(lon2d[m].mean()), float(lat2d[m].mean())
        town = min(REFERENCE_TOWNS, key=lambda t: haversine_m(lon_c, lat_c, t[1], t[2]))[0]
        rows.append(dict(segment_id=sid, label=f"{sid:02d} {town}", lon_min=edges[b], lon_max=edges[b + 1],
                         lon_c=lon_c, lat_c=lat_c, n_cells=int(m.sum())))
    return CoastSegments(pd.DataFrame(rows), cell_segment, grid)
