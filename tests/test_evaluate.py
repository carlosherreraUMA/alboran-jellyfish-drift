import numpy as np
import pandas as pd
import pytest

from jellydrift.coast import build_segments, coastal_cells
from jellydrift.evaluate import daily_index, midrank_percentiles, reference_matrix, score_events
from jellydrift.field import Grid


def test_midrank_percentiles_handle_ties_and_nan():
    ref = np.array([[0.0, 0.0, 1.0, 2.0, np.nan]])
    pct = midrank_percentiles(ref)
    np.testing.assert_allclose(pct[0, :4], [0.25, 0.25, 0.625, 0.875])
    assert np.isnan(pct[0, 4])


def test_reference_matrix_lag_and_edges():
    dates = pd.date_range("2022-01-01", periods=10)
    index = pd.DataFrame({0: [0, 0, 5, 0, 0, 0, 0, 0, 0, 0], 1: 0.0}, index=dates, dtype=float)
    events = pd.DataFrame({"segment_id": [0, 0], "date": pd.to_datetime(["2022-01-04", "2022-01-01"])})
    ref = reference_matrix(index, events, window=2, lag=1)
    np.testing.assert_allclose(ref[0], [0, 5, 5, 0, 0])  # 5 on Jan 3 also counts on Jan 4 (lag 1)
    assert np.isnan(ref[1, :2]).all() and ref[1, 2] == 0


def test_score_is_high_for_perfect_model_and_half_for_random():
    rng = np.random.default_rng(0)
    n_events, window = 80, 15
    perfect = rng.random((n_events, 2 * window + 1))
    perfect[:, window] = 2.0  # event day always has the highest index
    res = score_events(perfect, n_permutations=2000)
    assert res["mean_score"] > 0.98 and res["p_value"] < 0.001

    noise = rng.random((n_events, 2 * window + 1))
    res = score_events(noise, n_permutations=2000)
    assert abs(res["mean_score"] - 0.5) < 0.08 and res["p_value"] > 0.01


def test_all_zero_index_gives_neutral_score():
    res = score_events(np.zeros((30, 31)), n_permutations=500)
    assert res["mean_score"] == pytest.approx(0.5)
    assert res["hit_rate"] == 0.0


def _toy_segments():
    lon = np.linspace(-5.0, -4.0, 21)
    lat = np.linspace(36.0, 37.0, 21)
    land = np.zeros((21, 21), dtype=bool)
    land[15:, :] = True  # coast running east-west at ~36.7N
    return build_segments(Grid(lon, lat), land, lat_min=36.0, lon_min=-5.0, lon_max=-4.0, width_deg=0.25)


def test_coastal_cells_and_segments():
    segs = _toy_segments()
    assert coastal_cells(np.pad(np.ones((1, 1), bool), 1)).sum() == 8
    assert len(segs.table) == 4
    assert (segs.cell_segment >= 0).sum() == 20  # one row of coastal cells; bins are [lon_min, lon_max)
    seg, dist = segs.nearest_segment(np.array([-4.9, -4.1]), np.array([36.66, 36.0]), max_distance_km=5)
    assert seg[0] == 0 and seg[1] == -1 and dist[1] > 50


def test_daily_index_normalisation():
    segs = _toy_segments()
    i, j = np.nonzero(segs.cell_segment == 1)
    st = pd.DataFrame({
        "source": pd.Categorical(["a", "a", "b"], categories=["a", "b"]),
        "strand_time": pd.to_datetime(["2022-01-02 03:00", "2022-01-02 20:00", "2022-01-03 01:00"]),
        "cell_i": [i[0]] * 3, "cell_j": [j[0]] * 3, "age_days": [5.0, 50.0, 1.0],
    })
    dates = pd.date_range("2022-01-01", periods=4)
    idx = daily_index(st, segs, dates, {"a": 100.0, "b": 400.0}, max_age_days=10)
    assert idx["a"].loc["2022-01-02", 1] == pytest.approx(10.0)  # 1 stranding (age<=10) per 100/day -> 10 per 1000
    assert idx["b"].loc["2022-01-03", 1] == pytest.approx(2.5)
    assert idx["all_sources"].loc["2022-01-02", 1] == pytest.approx(2.0)  # 1 per 500/day
    assert idx["a"].to_numpy().sum() == pytest.approx(10.0)
