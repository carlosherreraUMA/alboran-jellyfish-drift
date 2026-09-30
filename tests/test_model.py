import numpy as np
import pytest

from jellydrift.field import Grid, GriddedLayer, LayeredForcing
from jellydrift.geo import EARTH_RADIUS_M, haversine_m, solar_elevation_deg
from jellydrift.model import DriftParams, Simulation, Source

T0 = np.datetime64("2022-03-01T00:00:00", "s")
LON = np.linspace(-5.0, -3.0, 41)  # 0.05 deg; linspace keeps exact values like -3.5
LAT = np.linspace(35.0, 37.0, 41)


def hourly_times(days):
    return T0 + np.arange(int(days * 24) + 1) * np.timedelta64(3600, "s")


def constant_layer(u0, v0, days=10, land=None, name="surface", depth=0.0):
    times = hourly_times(days)
    u = np.full((len(times), len(LAT), len(LON)), u0, dtype=np.float32)
    v = np.full_like(u, v0)
    if land is not None:
        u[:, land] = np.nan
        v[:, land] = np.nan
    return GriddedLayer.from_arrays(times, LON, LAT, u, v, name=name, depth=depth)


def test_grid_bilinear_is_exact_for_linear_fields():
    grid = Grid(LON, LAT)
    lon2d, lat2d = np.meshgrid(LON, LAT)
    f = 3.0 * lon2d - 2.0 * lat2d + 1.0
    pts_lon = np.array([-4.123, -3.51, -4.999])
    pts_lat = np.array([35.77, 36.02, 36.99])
    np.testing.assert_allclose(grid.bilinear(f, pts_lon, pts_lat), 3 * pts_lon - 2 * pts_lat + 1, atol=1e-9)


def test_irregular_grid_rejected():
    with pytest.raises(ValueError):
        Grid(np.array([0.0, 1.0, 3.0]), LAT)


def test_linear_time_interpolation():
    times = T0 + np.array([0, 3600], dtype="timedelta64[s]")
    u = np.stack([np.zeros((len(LAT), len(LON))), np.ones((len(LAT), len(LON)))]).astype(np.float32)
    layer = GriddedLayer.from_arrays(times, LON, LAT, u, np.zeros_like(u))
    uu, _ = layer.velocity(T0 + np.timedelta64(900, "s"), np.array([-4.0]), np.array([36.0]))
    assert uu[0] == pytest.approx(0.25)
    with pytest.raises(ValueError):
        layer.velocity(T0 + np.timedelta64(7200, "s"), np.array([-4.0]), np.array([36.0]))


def test_uniform_flow_displacement_matches_physics():
    sim = Simulation(LayeredForcing([constant_layer(0.5, 0.2)]), [], DriftParams(dt_seconds=900))
    lons, lats = sim.integrate([-4.5], [36.0], T0, T0 + np.timedelta64(86400, "s"))
    east = haversine_m(-4.5, lats[-1, 0], lons[-1, 0], lats[-1, 0])
    north = (lats[-1, 0] - 36.0) * np.pi / 180 * EARTH_RADIUS_M
    assert north == pytest.approx(0.2 * 86400, rel=1e-6)
    assert east == pytest.approx(0.5 * 86400, rel=5e-3)  # meridian convergence over 17 km of drift


def test_solid_body_rotation_closes_orbit():
    omega = 1.5e-5  # rad/s -> period ~4.8 days
    lon_c, lat_c = -4.0, 36.0
    lon2d, lat2d = np.meshgrid(LON, LAT)
    deg = np.pi / 180 * EARTH_RADIUS_M
    u2d = -omega * (lat2d - lat_c) * deg * np.cos(np.deg2rad(lat2d))
    v2d = omega * (lon2d - lon_c) * deg
    times = hourly_times(6)
    u = np.broadcast_to(u2d, (len(times), *u2d.shape)).astype(np.float32)
    v = np.broadcast_to(v2d, (len(times), *v2d.shape)).astype(np.float32)
    layer = GriddedLayer.from_arrays(times, LON, LAT, u, v)
    period = int(round(2 * np.pi / omega / 900)) * 900
    sim = Simulation(LayeredForcing([layer]), [], DriftParams(dt_seconds=900))
    lons, lats = sim.integrate([-4.3], [36.0], T0, T0 + np.timedelta64(period, "s"))
    radius = np.hypot(lons[:, 0] - lon_c, lats[:, 0] - lat_c)
    assert np.ptp(radius) < 1e-3  # degrees
    assert haversine_m(lons[-1, 0], lats[-1, 0], -4.3, 36.0) < 2_000  # < 2 km after 33 km-radius orbit


def test_stranding_time_and_location_on_eastern_wall():
    land = np.broadcast_to(LON >= -3.5, (len(LAT), len(LON)))
    forcing = LayeredForcing([constant_layer(1.0, 0.0, land=land)])
    src = Source("point", -4.5, -4.5, 36.0, 36.0, particles_per_day=96)  # one particle per 15-min step
    sim = Simulation(forcing, [src], DriftParams(dt_seconds=900, horizontal_diffusivity=0.0))
    res = sim.run(T0, T0 + np.timedelta64(3 * 86400, "s"), release_until=T0 + np.timedelta64(900, "s"))
    assert len(res.strandings) == 1
    row = res.strandings.iloc[0]
    # 0.975 deg of longitude at 36N (to the land-cell edge) is 87.7 km -> ~24.4 h at 1 m/s,
    # plus a short slow-down in the last half cell where velocity is interpolated towards 0.
    hours = (row["strand_time"] - T0) / np.timedelta64(1, "h")
    assert 24.0 < hours < 26.0
    assert LON[row["cell_j"]] == pytest.approx(-3.55)
    s = res.summary.iloc[0]
    assert s.released == s.stranded + s.exited + s.expired + s.active_at_end == 1


def test_particles_leaving_domain_are_counted_as_exited():
    forcing = LayeredForcing([constant_layer(-1.0, 0.0)])
    src = Source("west", -4.95, -4.9, 36.0, 36.1, particles_per_day=200)
    res = Simulation(forcing, [src], DriftParams(horizontal_diffusivity=0.0)).run(
        T0, T0 + np.timedelta64(2 * 86400, "s"), release_until=T0 + np.timedelta64(86400, "s"))
    s = res.summary.iloc[0]
    assert s.released == 200
    assert s.exited == 200 and s.stranded == 0


def test_release_avoids_land_and_runs_are_reproducible():
    land = np.broadcast_to((LAT >= 36.5)[:, None], (len(LAT), len(LON)))
    forcing = LayeredForcing([constant_layer(0.0, 0.3, land=land)])
    src = Source("box", -4.5, -3.5, 36.0, 37.0, particles_per_day=500)  # half the box is land
    runs = [Simulation(forcing, [src], DriftParams(seed=42, horizontal_diffusivity=5.0)).run(
        T0, T0 + np.timedelta64(3 * 86400, "s"), release_until=T0 + np.timedelta64(86400, "s")) for _ in range(2)]
    a, b = (r.strandings.sort_values("pid").reset_index(drop=True) for r in runs)
    assert len(a) > 0
    assert a.equals(b)
    assert (a["lat"] < 36.5).all()


def test_solar_elevation():
    noon = solar_elevation_deg(np.datetime64("2022-06-21T12:19:00"), -4.42, 36.72)
    assert noon == pytest.approx(90 - 36.72 + 23.44, abs=1.0)
    assert solar_elevation_deg(np.datetime64("2022-06-21T00:19:00"), -4.42, 36.72) < -20
