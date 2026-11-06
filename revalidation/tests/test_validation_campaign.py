import numpy as np
import pandas as pd
import pytest
from astropy.wcs import WCS

from revalidation.run_validation_campaign import (
    _aperture_metrics,
    _gaia_query,
    _ztf_summary_table,
    parse_delimited_response,
)


def test_parse_delimited_response_skips_service_comments():
    parsed = parse_delimited_response(
        "# service metadata\n# another comment\nid,value\n1,2.5\n",
        sep=",",
        expected_columns=("id", "value"),
    )
    assert len(parsed) == 1
    assert int(parsed.loc[0, "id"]) == 1
    assert float(parsed.loc[0, "value"]) == pytest.approx(2.5)


def test_gaia_query_contains_both_epoch_and_variability_tables():
    query, basic = _gaia_query(179.9, -0.05, 15.0 / 3600.0)
    assert "gaiadr3.gaia_source" in query
    assert "gaiadr3.vari_classifier_result" in query
    assert "gaiadr3.vari_summary" in query
    assert "CIRCLE('ICRS'" in query
    assert "gaiadr3.vari_classifier_result" not in basic


def test_ztf_summary_separates_strict_catflags_zero_photometry():
    candidate = {"ra_deg": 179.9, "dec_deg": -0.05}
    frame = pd.DataFrame({
        "oid": [123, 123, 123],
        "ra": [179.90001, 179.90001, 179.90001],
        "dec": [-0.05001, -0.05001, -0.05001],
        "mjd": [59000.0, 59001.0, 59002.0],
        "mag": [16.1, 16.2, 16.4],
        "magerr": [0.03, 0.03, 0.04],
        "filtercode": ["zg", "zg", "zr"],
        "catflags": [0, 32768, 0],
    })
    summary, result = _ztf_summary_table("1978", candidate, frame, "ztf_dr24")
    assert result["status"] == "COMPLETED"
    assert len(summary) == 1
    assert int(summary.loc[0, "rows_returned"]) == 3
    assert int(summary.loc[0, "strict_catflags_zero_rows"]) == 2
    assert np.isfinite(summary.loc[0, "strict_mag_reduced_chi2"])
    assert float(summary.loc[0, "nearest_separation_arcsec"]) < 1.0


def test_aperture_remeasurement_returns_mask_sensitivity_and_mag():
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [51.0, 51.0]
    wcs.wcs.crval = [180.0, 0.0]
    wcs.wcs.cdelt = [-2.75 / 3600.0, 2.75 / 3600.0]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = wcs.to_header()
    header["MAGZP"] = 20.0
    yy, xx = np.indices((101, 101))
    intensity = 10.0 + 150.0 * np.exp(-((xx - 50.0) ** 2 + (yy - 50.0) ** 2) / 4.0)
    uncertainty = np.ones((101, 101), dtype=float)
    mask = np.zeros((101, 101), dtype=np.uint32)
    mask[50, 50] = 2 ** 30
    result = _aperture_metrics(
        intensity, uncertainty, mask, header, 180.0, 0.0, 13.0
    )
    assert result["aperture_pixels_clean"] < result["aperture_pixels_raw"]
    assert result["recommended_masked"]["aperture_flux_dn"] > 0
    assert np.isfinite(result["recommended_masked"]["aperture_mag_uncorrected"])
    assert result["recommended_masked"]["aperture_flux_dn"] != pytest.approx(
        result["raw"]["aperture_flux_dn"]
    )
