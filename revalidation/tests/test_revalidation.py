import math

import numpy as np
import pandas as pd
import pytest

from revalidation.reproduce_neowise import (
    CANDIDATES,
    angular_separation_arcsec,
    arithmetic_center_reduced_chi_square,
    build_groups,
    leave_one_out_chi_square,
    prepare_detections,
    parse_dec_degrees,
    parse_ra_degrees,
    reduced_chi_square,
    weighted_mean,
)


def test_weighted_mean_uses_inverse_variance_weights():
    assert weighted_mean([1.0, 2.0], [1.0, 2.0]) == pytest.approx(1.2)


def test_reduced_chi_square_uses_weighted_mean():
    assert reduced_chi_square([1.0, 3.0], [1.0, 1.0]) == pytest.approx(2.0)


def test_weighted_center_can_differ_from_arithmetic_center():
    weighted = reduced_chi_square([0.0, 2.0], [0.1, 1.0])
    arithmetic = arithmetic_center_reduced_chi_square([0.0, 2.0], [0.1, 1.0])
    assert weighted != pytest.approx(arithmetic)


def test_reduced_chi_square_returns_nan_for_one_point():
    assert math.isnan(reduced_chi_square([1.0], [0.2]))


def test_invalid_uncertainties_are_excluded():
    assert weighted_mean([1.0, 99.0, 3.0], [1.0, 0.0, 1.0]) == pytest.approx(2.0)


def test_leave_one_out_recomputes_every_omission():
    result = leave_one_out_chi_square([10.0, 10.1, 10.2], [0.1, 0.1, 0.1])
    assert result.shape == (3,)
    assert np.all(np.isfinite(result))


def test_spherical_separation_handles_ra_wrap():
    separation = angular_separation_arcsec(359.999, 0.0, [0.001], [0.0])
    assert separation[0] == pytest.approx(7.2, rel=1e-6)


def _synthetic_detections() -> pd.DataFrame:
    rows = []
    for candidate_id, coordinate in CANDIDATES.items():
        for index in range(12):
            rows.append({
                "ra": coordinate["ra_deg"] + (index - 5.5) * 0.000005,
                "dec": coordinate["dec_deg"] + (index % 3 - 1) * 0.000004,
                "mjd": 59000.0 + index * 180.0,
                "w1mpro": 12.0 + (index % 4) * 0.02,
                "w1sigmpro": 0.03,
                "w1snr": 20.0,
                "w2mpro": 11.8 + (index % 4) * 0.02,
                "w2sigmpro": 0.04,
                "w2snr": 15.0,
                "qual_frame": 10,
                "ph_qual": "AAxx",
                "cc_flags": "0000",
                "scan_id": f"{44000 + index:05d}a",
                "frame_num": index + 1,
                "qi_fact": 1.0,
                "saa_sep": 10.0,
                "moon_masked": "00",
            })
    # This row lies at the previously excluded EXP-001 Candidate #1 position.
    rows.append({
        "ra": 179.936465,
        "dec": -0.0897475,
        "mjd": 59000.0,
        "w1mpro": 12.0,
        "w1sigmpro": 0.03,
        "w1snr": 20.0,
        "w2mpro": 11.8,
        "w2sigmpro": 0.04,
        "w2snr": 15.0,
        "qual_frame": 10,
        "ph_qual": "AAxx",
        "cc_flags": "0000",
        "scan_id": "44000a",
        "frame_num": 1,
        "qi_fact": 1.0,
        "saa_sep": 10.0,
        "moon_masked": "00",
    })
    return pd.DataFrame(rows)


def test_prepare_detections_preserves_original_and_expanded_quality_views():
    prepared, counts = prepare_detections(_synthetic_detections())
    assert counts["raw_rows"] == 25
    assert counts["numerically_usable_rows"] == 25
    assert counts["candidate_1_excluded_rows"] == 1
    assert len(prepared) == 24
    assert prepared["original_good_quality"].all()
    assert prepared["expanded_quality_available_filters"].all()
    assert counts["expanded_quality_filter_complete"] is True


def test_coordinate_candidate_matching_survives_group_label_changes():
    prepared, _ = prepare_detections(_synthetic_detections())
    _, groups = build_groups(prepared, 1.5)
    assert len(groups) == 2
    for candidate_id, coordinate in CANDIDATES.items():
        separation = angular_separation_arcsec(
            coordinate["ra_deg"], coordinate["dec_deg"],
            groups["ra"].to_numpy(), groups["dec"].to_numpy(),
        )
        assert float(np.min(separation)) < 1.0


def test_parse_ra_sexagesimal_hours():
    assert parse_ra_degrees("12 00 00.0") == pytest.approx(180.0)


def test_parse_ra_decimal_degrees():
    assert parse_ra_degrees("180.25") == pytest.approx(180.25)


def test_parse_negative_sexagesimal_declination():
    assert parse_dec_degrees("-00 05 24.0") == pytest.approx(-0.09)


def test_parse_negative_decimal_declination_preserves_sign():
    assert parse_dec_degrees("-0.09") == pytest.approx(-0.09)
