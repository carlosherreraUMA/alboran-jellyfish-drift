"""Geographic helpers: metric <-> degree conversions, distances and solar elevation."""

from __future__ import annotations

import numpy as np

EARTH_RADIUS_M = 6_371_000.0
DEG_PER_M_LAT = 180.0 / (np.pi * EARTH_RADIUS_M)


def meters_to_degrees(dx_m, dy_m, lat_deg):
    """Convert eastward/northward displacements in metres to (dlon, dlat) in degrees."""
    dlat = dy_m * DEG_PER_M_LAT
    dlon = dx_m * DEG_PER_M_LAT / np.cos(np.deg2rad(lat_deg))
    return dlon, dlat


def velocity_to_degrees_per_second(u_ms, v_ms, lat_deg):
    """Convert (u, v) in m/s to (dlon/dt, dlat/dt) in degrees per second."""
    return meters_to_degrees(u_ms, v_ms, lat_deg)


def haversine_m(lon1, lat1, lon2, lat2):
    """Great-circle distance in metres (broadcasts over arrays)."""
    lon1, lat1, lon2, lat2 = map(np.deg2rad, (lon1, lat1, lon2, lat2))
    a = (
        np.sin((lat2 - lat1) / 2) ** 2
        + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    )
    return 2 * EARTH_RADIUS_M * np.arcsin(np.sqrt(a))


_J2000 = np.datetime64("2000-01-01T12:00:00", "s")


def solar_elevation_deg(time_utc, lon_deg, lat_deg):
    """Approximate solar elevation angle (degrees) for UTC time(s).

    Uses the low-precision US Naval Observatory algorithm (accurate to ~1 degree),
    which is more than enough to separate day from night for diel vertical migration.
    """
    t = np.asarray(time_utc, dtype="datetime64[s]")
    d = (t - _J2000).astype("float64") / 86400.0
    g = np.deg2rad((357.529 + 0.98560028 * d) % 360.0)
    q = (280.459 + 0.98564736 * d) % 360.0
    ecl_lon = np.deg2rad(q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    obliq = np.deg2rad(23.439 - 0.00000036 * d)
    ra = np.arctan2(np.cos(obliq) * np.sin(ecl_lon), np.cos(ecl_lon))
    decl = np.arcsin(np.sin(obliq) * np.sin(ecl_lon))
    gmst_h = (18.697374558 + 24.06570982441908 * d) % 24.0
    hour_angle = np.deg2rad(gmst_h * 15.0 + lon_deg) - ra
    lat = np.deg2rad(lat_deg)
    sin_elev = np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(hour_angle)
    return np.rad2deg(np.arcsin(np.clip(sin_elev, -1.0, 1.0)))
