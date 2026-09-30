import numpy as np
import pytest

from jellydrift.behaviour import DielMigration, FixedDepth, behaviour_from_config
from jellydrift.field import Bathymetry, Grid, GriddedLayer, LayeredForcing
from jellydrift.model import DriftParams, Simulation, Source

T0 = np.datetime64("2022-03-01T00:00:00", "s")
LON = np.linspace(-5.0, -3.0, 41)
LAT = np.linspace(35.0, 37.0, 41)
NOON = np.datetime64("2022-03-03T12:15:00", "s")  # solar noon at 4W is ~12:28 UTC in March
MIDNIGHT = np.datetime64("2022-03-03T00:15:00", "s")


def layer(u0, depth, land=None, days=10):
    times = T0 + np.arange(days * 24 + 1) * np.timedelta64(3600, "s")
    u = np.full((len(times), len(LAT), len(LON)), u0, dtype=np.float32)
    if land is not None:
        u[:, land] = np.nan
    return GriddedLayer.from_arrays(times, LON, LAT, u, np.where(np.isnan(u), np.nan, 0.0), depth=depth,
                                    name=f"z{depth}")


def test_vertical_interpolation_between_levels_and_below_seabed():
    shelf = np.zeros((len(LAT), len(LON)), dtype=bool)
    shelf[:, LON > -3.6] = True  # the 100 m level does not exist east of 3.6W
    f = LayeredForcing([layer(-1.0, 100.0, land=shelf), layer(1.0, 0.0)])  # order does not matter
    assert list(f.depths) == [0.0, 100.0]
    lon = np.array([-4.5, -4.5, -4.5, -3.2])
    lat = np.full(4, 36.0)
    u, _ = f.velocity(NOON, lon, lat, np.array([0.0, 25.0, 150.0, 50.0]))
    np.testing.assert_allclose(u, [1.0, 0.5, -1.0, 1.0])  # linear in depth; deepest level below; shelf -> upper
    np.testing.assert_allclose(f.velocity(NOON, lon, lat)[0], 1.0)  # no depth -> surface


def test_bottom_depth_from_bathymetry_and_proxy():
    shelf = np.zeros((len(LAT), len(LON)), dtype=bool)
    shelf[:, LON > -3.6] = True
    f = LayeredForcing([layer(0.0, 0.0), layer(0.0, 150.0, land=shelf)])
    np.testing.assert_allclose(f.bottom_depth(np.array([-4.5, -3.2]), np.array([36.0, 36.0])), [np.inf, 150.0])
    depth = np.full((len(LAT), len(LON)), 800.0)
    f2 = LayeredForcing(f.layers, Bathymetry(Grid(LON, LAT), depth))
    np.testing.assert_allclose(f2.bottom_depth(np.array([-4.5]), np.array([36.0])), [800.0])


def test_diel_migration_targets_follow_sun_and_seabed():
    b = DielMigration(day_depth_m=150, night_depth_m=5, bottom_clearance_m=100, no_migration_below_m=100)
    lon = np.full(4, -4.0)
    lat = np.full(4, 36.0)
    bottom = np.array([2000.0, 200.0, 120.0, 80.0])
    day = b.target_depth(NOON, lon, lat, np.full(4, 150.0), bottom)
    np.testing.assert_allclose(day, [150.0, 100.0, 20.0, 5.0])  # open sea, 100 m above bottom, shallow: none
    night = b.target_depth(MIDNIGHT, lon, lat, np.full(4, 150.0), bottom)
    np.testing.assert_allclose(night, 5.0)


def test_swimming_speed_limits_vertical_motion():
    b = DielMigration(speed_m_per_h=100)
    z = b.swim(np.array([5.0, 150.0, 60.0]), np.array([150.0, 5.0, 70.0]), dt_seconds=1800)
    np.testing.assert_allclose(z, [55.0, 100.0, 70.0])


def test_individual_day_depths_variability_and_reproducibility():
    b = DielMigration(day_depth_m=150, day_depth_rel_sd=0.2)
    d1 = b.individual_day_depths(20000, np.random.default_rng(3))
    d2 = b.individual_day_depths(20000, np.random.default_rng(3))
    np.testing.assert_array_equal(d1, d2)
    assert d1.mean() == pytest.approx(150, rel=0.01) and d1.std() == pytest.approx(30, rel=0.05)
    assert d1.min() >= 5 and d1.max() <= 300
    np.testing.assert_array_equal(DielMigration(day_depth_m=150).individual_day_depths(3, None), 150.0)


def test_behaviour_from_config():
    assert isinstance(behaviour_from_config(None), FixedDepth)
    b = behaviour_from_config(300, {"speed_m_per_h": 50})
    assert b.day_depth_m == 300 and b.speed_m_per_h == 50


def test_migrating_particles_feel_deep_currents_by_day():
    """Surface flow +0.5 m/s east, 150 m flow -0.5 m/s: migrants drift less than surface drifters."""
    f = LayeredForcing([layer(0.5, 0.0), layer(-0.5, 150.0)], Bathymetry(Grid(LON, LAT), np.full((41, 41), 2000.0)))
    src = Source("point", -4.5, -4.5, 36.0, 36.0, particles_per_day=48)  # one particle per 30-min step
    kw = dict(start=T0, end=T0 + np.timedelta64(2 * 86400, "s"), release_until=T0 + np.timedelta64(1800, "s"),
              snapshot_every_hours=0.5)
    params = DriftParams(dt_seconds=1800, horizontal_diffusivity=0.0)
    passive = Simulation(f, [src], params).run(**kw).snapshots
    migr = Simulation(f, [src], params, DielMigration(day_depth_m=150, speed_m_per_h=100)).run(**kw).snapshots
    dx_passive = passive[-1].lon.iloc[0] - (-4.5)
    dx_migr = migr[-1].lon.iloc[0] - (-4.5)
    assert dx_passive > 0.8  # ~86 km east in 2 days
    assert dx_migr < 0.5 * dx_passive
    depths = np.array([s.z.iloc[0] for s in migr])
    assert depths.min() == pytest.approx(5.0) and depths.max() == pytest.approx(150.0)
    assert np.all(np.abs(np.diff(depths)) <= 50.0 + 1e-9)  # 100 m/h over 30-min snapshots


def test_deep_animals_do_not_strand_on_steep_coasts():
    wall = np.broadcast_to(LON >= -3.5, (41, 41))
    f = LayeredForcing([layer(1.0, 0.0, land=wall)])
    src = Source("point", -3.7, -3.7, 36.0, 36.0, particles_per_day=48)
    kw = dict(start=T0, end=T0 + np.timedelta64(86400, "s"), release_until=T0 + np.timedelta64(1800, "s"))
    params = DriftParams(dt_seconds=1800, horizontal_diffusivity=0.0, max_stranding_depth_m=25)
    shallow = Simulation(f, [src], params, FixedDepth(5.0)).run(**kw)
    deep = Simulation(f, [src], params, FixedDepth(150.0)).run(**kw)
    assert shallow.summary.stranded.iloc[0] == 1
    assert deep.summary.stranded.iloc[0] == 0 and deep.summary.active_at_end.iloc[0] == 1


def test_stacked_and_per_level_interpolation_agree():
    rng = np.random.default_rng(5)
    times = T0 + np.arange(49) * np.timedelta64(3600, "s")
    shelf = np.zeros((41, 41), dtype=bool)
    shelf[:, LON > -3.4] = True
    layers = []
    for depth in (0.0, 10.0, 50.0, 150.0):
        u = rng.normal(0, 0.3, (49, 41, 41)).astype(np.float32)
        v = rng.normal(0, 0.3, (49, 41, 41)).astype(np.float32)
        if depth >= 50:
            u[:, shelf] = np.nan
            v[:, shelf] = np.nan
        layers.append(GriddedLayer.from_arrays(times, LON, LAT, u, v, depth=depth))
    fast = LayeredForcing(layers)
    slow = LayeredForcing(layers)
    slow.stacked = False
    assert fast.stacked
    n = 5000
    lon, lat, z = rng.uniform(-4.99, -3.01, n), rng.uniform(35.01, 36.99, n), rng.uniform(0, 200, n)
    t = T0 + np.timedelta64(5 * 3600 + 900, "s")
    for a, b in zip(fast.velocity(t, lon, lat, z), slow.velocity(t, lon, lat, z), strict=True):
        np.testing.assert_allclose(a, b, atol=1e-6)
