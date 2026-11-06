"""Independent validation campaign for NEOWISE candidate revalidation v2.

Tracks:
1. GCVS positive-control cross-match (VizieR mirrors, spatially constrained query).
2. Gaia DR3 astrometry/variability and ZTF public optical light curves.
3. Direct single-exposure NEOWISE FITS cutouts (intensity/uncertainty/mask)
   and repeatable aperture photometry on the eight most influential frames.
4. Cadence-preserving synthetic sinusoid injection/recovery on internal stable
   controls, plus residual-bootstrap null experiments.
5. A spatially held-out NEOWISE field with thresholds frozen from EXP-005.

None of these checks is, by itself, a discovery claim. In particular, the
held-out field is independent sky area but still the same NEOWISE instrument
and archive, and the internal stable controls are not externally proven
non-variables.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import re
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.wcs import WCS

try:
    from .reproduce_neowise import (
        CANDIDATES,
        DEC_MAX,
        DEC_MIN,
        GROUP_RADIUS_ARCSEC,
        MIN_OBSERVATIONS,
        RA_MAX,
        RA_MIN,
        angular_separation_arcsec,
        build_groups,
        parse_dec_degrees,
        parse_ra_degrees,
        prepare_detections,
        reduced_chi_square,
        retrieve_region,
        weighted_mean,
    )
except ImportError:  # Direct execution as "python revalidation/run_validation_campaign.py"
    from reproduce_neowise import (
        CANDIDATES,
        DEC_MAX,
        DEC_MIN,
        GROUP_RADIUS_ARCSEC,
        MIN_OBSERVATIONS,
        RA_MAX,
        RA_MIN,
        angular_separation_arcsec,
        build_groups,
        parse_dec_degrees,
        parse_ra_degrees,
        prepare_detections,
        reduced_chi_square,
        retrieve_region,
        weighted_mean,
    )

GCVS_TABLE = "B/gcvs/gcvs_cat"
GCVS_MIRRORS = (
    "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv",
    "https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
    "https://vizier.nao.ac.jp/viz-bin/asu-tsv",
)
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
ZTF_LIGHTCURVE_API = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
NEOWISE_IMAGE_ROOT = "https://irsa.ipac.caltech.edu/ibe/data/wise/neowiser/p1bm_frm"

GAIA_RADIUS_DEG = 15.0 / 3600.0
ZTF_RADIUS_DEG = 5.0 / 3600.0
IMAGE_PRODUCTS = ("int", "unc", "msk")
IMAGE_BANDS = (1, 2)
IMAGE_TOP_N_PER_CANDIDATE = 4
APERTURE_RADIUS_PIX = 3.0
ANNULUS_INNER_PIX = 6.0
ANNULUS_OUTER_PIX = 12.0
# Conservative set of NEOWISE mask bits recommended for coaddition. We report
# both all finite pixels and this conservative-mask sensitivity measurement.
NEOWISE_RECOMMENDED_MASK_BITS = 1076887130

INJECTION_AMPLITUDES_MAG = (0.02, 0.05, 0.10, 0.20)
INJECTION_PERIODS_DAYS = (30.0, 60.0, 120.0, 180.0, 365.0)
INJECTION_PHASES = 12
NULL_BOOTSTRAP_REPLICATES = 500
VARIABILITY_THRESHOLD = 2.0
RANDOM_SEED = 20261106

# A separate sky field, 2 degrees from the development field in RA.
# It does not overlap the EXP-005 field.
HELDOUT_BOUNDS = (182.0, 182.2, -0.1, 0.1)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if value is pd.NA:
        return None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(value: dict[str, Any], path: Path) -> None:
    path.write_text(
        json.dumps(json_safe(value), indent=2, allow_nan=False),
        encoding="utf-8",
    )


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.12g")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_delimited_response(
    text: str,
    *,
    sep: str = ",",
    expected_columns: tuple[str, ...] = (),
) -> pd.DataFrame:
    """Parse a text table while ignoring service comment/metadata lines."""
    lines = [
        line for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        return pd.DataFrame(columns=list(expected_columns))
    frame = pd.read_csv(io.StringIO("\n".join(lines)), sep=sep)
    frame.columns = [str(c).strip().lstrip("\ufeff") for c in frame.columns]
    if expected_columns and not set(expected_columns).issubset(frame.columns):
        raise ValueError(
            f"Response lacks expected columns {expected_columns}; got {list(frame.columns)}"
        )
    return frame


def _candidate_frame(
    candidate_id: str,
    candidate: dict[str, float],
    detections: pd.DataFrame,
    groups: pd.DataFrame,
) -> pd.DataFrame:
    raw_sep = angular_separation_arcsec(
        candidate["ra_deg"], candidate["dec_deg"],
        detections["ra"].to_numpy(dtype=float),
        detections["dec"].to_numpy(dtype=float),
    )
    source_sep = angular_separation_arcsec(
        candidate["ra_deg"], candidate["dec_deg"],
        groups["ra"].to_numpy(dtype=float),
        groups["dec"].to_numpy(dtype=float),
    )
    closest_source_idx = int(np.argmin(source_sep))
    closest_source = groups.iloc[closest_source_idx]
    closest_detection_idx = int(np.argmin(raw_sep))
    closest_detection = detections.iloc[closest_detection_idx]
    return pd.DataFrame([{
        "candidate_id": candidate_id,
        "candidate_ra_deg": candidate["ra_deg"],
        "candidate_dec_deg": candidate["dec_deg"],
        "nearest_neowise_detection_sep_arcsec": float(raw_sep[closest_detection_idx]),
        "nearest_source_group_key": str(closest_source["source_key"]),
        "nearest_source_group_sep_arcsec": float(source_sep[closest_source_idx]),
        "source_group_match_5arcsec": bool(source_sep[closest_source_idx] <= 5.0),
        "source_group_n_observations": int(closest_source["n_observations"]),
        "source_group_w1_reduced_chi2": float(closest_source["w1_reduced_chi2"]),
        "nearest_detection_mjd": float(closest_detection["mjd"]),
    }])


def run_gcvs_check(
    detections: pd.DataFrame,
    groups: pd.DataFrame,
    out_dir: Path,
) -> dict[str, Any]:
    """Query GCVS using multiple VizieR mirrors and cross-match by sky position."""
    reference_path = out_dir / "gcvs_reference_population.csv"
    matches_path = out_dir / "gcvs_neowise_crossmatch.csv"
    errors: list[str] = []
    catalogue: pd.DataFrame | None = None
    successful_url: str | None = None

    spatial_params = {
        "-source": GCVS_TABLE,
        "-out": "GCVS,RAJ2000,DEJ2000,VarType",
        "-c": "180 0",
        "-c.rs": "12",
        "-out.max": "2000",
    }
    # Cone syntax varies between VizieR mirrors. If a valid cone query returns
    # no rows, retrieve the full catalogue (about 60k rows) and spatially filter
    # it locally, matching the method used in the original EXP-006 notebook.
    full_params = {
        "-source": GCVS_TABLE,
        "-out": "GCVS,RAJ2000,DEJ2000,VarType",
        "-out.max": "100000",
    }
    search_mode = "spatial_cone"
    for endpoint in GCVS_MIRRORS:
        try:
            response = requests.get(
                endpoint,
                params=spatial_params,
                timeout=(15, 55),
                headers={"User-Agent": "NEOWISE-revalidation/2.0 research"},
            )
            response.raise_for_status()
            parsed = parse_delimited_response(
                response.text,
                sep="\\t",
                expected_columns=("GCVS", "RAJ2000", "DEJ2000", "VarType"),
            )
            if parsed.empty:
                response = requests.get(
                    endpoint,
                    params=full_params,
                    timeout=(15, 90),
                    headers={"User-Agent": "NEOWISE-revalidation/2.0 research"},
                )
                response.raise_for_status()
                parsed = parse_delimited_response(
                    response.text,
                    sep="\\t",
                    expected_columns=("GCVS", "RAJ2000", "DEJ2000", "VarType"),
                )
                search_mode = "full_catalogue_local_spatial_filter"
            if parsed.empty:
                errors.append(
                    f"{urlparse(endpoint).netloc}: spatial and full catalogue queries returned no rows"
                )
                continue
            catalogue = parsed
            successful_url = response.url
            break
        except Exception as exc:
            errors.append(f"{urlparse(endpoint).netloc}: {type(exc).__name__}: {exc}")

    if catalogue is None:
        write_csv(
            pd.DataFrame(columns=[
                "object_id", "ra_deg", "dec_deg", "variability_type",
            ]),
            reference_path,
        )
        write_csv(
            pd.DataFrame(columns=[
                "object_id", "variability_type", "nearest_detection_sep_arcsec",
                "nearest_source_group_sep_arcsec", "source_group_key",
                "source_group_n_observations", "source_group_w1_reduced_chi2",
            ]),
            matches_path,
        )
        return {
            "status": "UNAVAILABLE",
            "catalogue": "GCVS 5.1",
            "errors": errors,
            "historical_notebook_result": {
                "catalogue_objects_in_field": 1,
                "matched_neowise_objects": 1,
                "variability_type": "RRAB",
                "nearest_detection_separation_arcsec": 1.2322,
                "validator_result": "VARIABLE_HIGH_QUALITY",
                "source_experiment": "EXP-006 / EXP-007 outputs in public notebook",
            },
            "interpretation": (
                "Current VizieR mirrors did not respond. The historical public "
                "notebook records one RRAB positive control and its validator "
                "classification; that historical result is preserved but is "
                "not counted as a fresh cross-match in this run."
            ),
            "reference_file": reference_path.name,
            "matches_file": matches_path.name,
        }

    # Force numeric dtypes explicitly: pandas can preserve object dtype when
    # the VizieR response contains masked coordinate cells or mixed tokens.
    catalogue["ra_deg"] = pd.to_numeric(
        catalogue["RAJ2000"].apply(parse_ra_degrees), errors="coerce"
    ).astype("float64")
    catalogue["dec_deg"] = pd.to_numeric(
        catalogue["DEJ2000"].apply(parse_dec_degrees), errors="coerce"
    ).astype("float64")
    valid_coordinate_mask = (
        np.isfinite(catalogue["ra_deg"].to_numpy(dtype=np.float64))
        & np.isfinite(catalogue["dec_deg"].to_numpy(dtype=np.float64))
    )
    valid = catalogue.loc[valid_coordinate_mask].copy()
    in_field = valid[
        valid["ra_deg"].between(RA_MIN, RA_MAX)
        & valid["dec_deg"].between(DEC_MIN, DEC_MAX)
    ].copy()
    reference = (
        in_field[["GCVS", "ra_deg", "dec_deg", "VarType"]]
        .drop_duplicates("GCVS")
        .rename(columns={
            "GCVS": "object_id",
            "VarType": "variability_type",
        })
        .reset_index(drop=True)
    )
    write_csv(reference, reference_path)

    rows: list[dict[str, Any]] = []
    for _, obj in reference.iterrows():
        dsep = angular_separation_arcsec(
            float(obj["ra_deg"]), float(obj["dec_deg"]),
            detections["ra"].to_numpy(dtype=float),
            detections["dec"].to_numpy(dtype=float),
        )
        gsep = angular_separation_arcsec(
            float(obj["ra_deg"]), float(obj["dec_deg"]),
            groups["ra"].to_numpy(dtype=float),
            groups["dec"].to_numpy(dtype=float),
        )
        di = int(np.argmin(dsep))
        gi = int(np.argmin(gsep))
        group = groups.iloc[gi]
        rows.append({
            "object_id": str(obj["object_id"]),
            "variability_type": str(obj["variability_type"]),
            "reference_ra_deg": float(obj["ra_deg"]),
            "reference_dec_deg": float(obj["dec_deg"]),
            "nearest_detection_sep_arcsec": float(dsep[di]),
            "nearest_detection_match_5arcsec": bool(dsep[di] <= 5.0),
            "nearest_source_group_sep_arcsec": float(gsep[gi]),
            "source_group_key": str(group["source_key"]),
            "source_group_match_5arcsec": bool(gsep[gi] <= 5.0),
            "source_group_n_observations": int(group["n_observations"]),
            "source_group_w1_reduced_chi2": float(group["w1_reduced_chi2"]),
        })
    matches = pd.DataFrame(rows)
    write_csv(matches, matches_path)

    return {
        "status": "COMPLETED",
        "catalogue": "GCVS 5.1",
        "service": successful_url,
        "search_mode": search_mode,
        "catalogue_rows_returned": int(len(catalogue)),
        "valid_coordinate_rows": int(len(valid)),
        "objects_inside_development_field": int(len(reference)),
        "objects_matched_to_neowise_within_5arcsec": (
            int(matches["nearest_detection_match_5arcsec"].sum())
            if len(matches) else 0
        ),
        "objects_matched_to_source_groups_within_5arcsec": (
            int(matches["source_group_match_5arcsec"].sum())
            if len(matches) else 0
        ),
        "interpretation": (
            "Positive-control consistency check only. The local GCVS sample is "
            "too small to estimate recall, precision, or false-positive rate."
        ),
        "reference_file": reference_path.name,
        "matches_file": matches_path.name,
    }


def _gaia_query(ra_deg: float, dec_deg: float, radius_deg: float) -> tuple[str, bool]:
    """Return ADQL and whether the query includes joined variability products."""
    rich = f"""
SELECT TOP 500
    g.source_id, g.ra, g.dec, g.pmra, g.pmdec, g.parallax,
    g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
    g.phot_variable_flag, g.ruwe,
    c.best_class_name, c.best_class_score,
    v.num_selected_g_fov, v.mean_obs_time_g_fov, v.time_duration_g_fov,
    v.mean_mag_g_fov, v.range_mag_g_fov
FROM gaiadr3.gaia_source AS g
LEFT OUTER JOIN gaiadr3.vari_classifier_result AS c
    ON g.source_id = c.source_id
LEFT OUTER JOIN gaiadr3.vari_summary AS v
    ON g.source_id = v.source_id
WHERE 1 = CONTAINS(
    POINT('ICRS', g.ra, g.dec),
    CIRCLE('ICRS', {ra_deg:.10f}, {dec_deg:.10f}, {radius_deg:.10f})
)
"""
    basic = f"""
SELECT TOP 500
    g.source_id, g.ra, g.dec, g.pmra, g.pmdec, g.parallax,
    g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
    g.phot_variable_flag, g.ruwe
FROM gaiadr3.gaia_source AS g
WHERE 1 = CONTAINS(
    POINT('ICRS', g.ra, g.dec),
    CIRCLE('ICRS', {ra_deg:.10f}, {dec_deg:.10f}, {radius_deg:.10f})
)
"""
    return rich, basic


def run_gaia_check(
    detections: pd.DataFrame,
    out_dir: Path,
) -> dict[str, Any]:
    """Cross-match Gaia DR3 astrometry and variability summary near candidates."""
    all_rows: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    errors: list[str] = []

    for candidate_id, candidate in CANDIDATES.items():
        candidate_detections = detections.loc[
            detections["legacy_candidate_id"].astype(str) == candidate_id
        ]
        if candidate_detections.empty:
            summary_rows.append({
                "candidate_id": candidate_id,
                "status": "NO_NEOWISE_DETECTION_ROWS",
            })
            continue

        median_mjd = float(candidate_detections["mjd"].median())
        wise_epoch_year = 2000.0 + (median_mjd - 51544.5) / 365.25
        rich_query, basic_query = _gaia_query(
            candidate["ra_deg"], candidate["dec_deg"], GAIA_RADIUS_DEG
        )
        payload = {
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "csv",
            "QUERY": rich_query,
        }
        result: pd.DataFrame | None = None
        mode = "variability_join"
        try:
            response = requests.post(
                GAIA_TAP,
                data=payload,
                timeout=(20, 150),
                headers={"User-Agent": "NEOWISE-revalidation/2.0 research"},
            )
            response.raise_for_status()
            result = parse_delimited_response(response.text, sep=",")
            if "source_id" not in result.columns:
                raise ValueError("Gaia TAP response did not contain source_id.")
        except Exception as exc:
            errors.append(
                f"Gaia {candidate_id} variability join: {type(exc).__name__}: {exc}"
            )
            mode = "basic_astrometry_fallback"
            try:
                response = requests.post(
                    GAIA_TAP,
                    data={**payload, "QUERY": basic_query},
                    timeout=(20, 120),
                    headers={"User-Agent": "NEOWISE-revalidation/2.0 research"},
                )
                response.raise_for_status()
                result = parse_delimited_response(response.text, sep=",")
                if "source_id" not in result.columns:
                    raise ValueError("Gaia fallback response did not contain source_id.")
            except Exception as fallback_exc:
                errors.append(
                    f"Gaia {candidate_id} basic query: "
                    f"{type(fallback_exc).__name__}: {fallback_exc}"
                )

        if result is None:
            summary_rows.append({
                "candidate_id": candidate_id,
                "status": "UNAVAILABLE",
                "query_mode": mode,
            })
            continue
        if result.empty:
            summary_rows.append({
                "candidate_id": candidate_id,
                "status": "NO_GAIA_SOURCES_WITHIN_15ARCSEC",
                "query_mode": mode,
            })
            continue

        result["candidate_id"] = candidate_id
        result["gaia_epoch_separation_arcsec"] = angular_separation_arcsec(
            candidate["ra_deg"], candidate["dec_deg"],
            pd.to_numeric(result["ra"], errors="coerce").to_numpy(dtype=float),
            pd.to_numeric(result["dec"], errors="coerce").to_numpy(dtype=float),
        )
        ra = pd.to_numeric(result["ra"], errors="coerce").to_numpy(dtype=float)
        dec = pd.to_numeric(result["dec"], errors="coerce").to_numpy(dtype=float)
        pmra = (
            pd.to_numeric(result["pmra"], errors="coerce").to_numpy(dtype=float)
            if "pmra" in result.columns else np.full(len(result), np.nan)
        )
        pmdec = (
            pd.to_numeric(result["pmdec"], errors="coerce").to_numpy(dtype=float)
            if "pmdec" in result.columns else np.full(len(result), np.nan)
        )
        dt = wise_epoch_year - 2016.0
        cosdec = np.cos(np.deg2rad(dec))
        with np.errstate(divide="ignore", invalid="ignore"):
            propagated_ra = (ra + (pmra * dt) / (3.6e6 * cosdec)) % 360.0
        propagated_dec = dec + (pmdec * dt) / 3.6e6
        valid_motion = (
            np.isfinite(propagated_ra) & np.isfinite(propagated_dec)
            & (np.abs(dec) < 89.9)
        )
        propagated_ra[~valid_motion] = ra[~valid_motion]
        propagated_dec[~valid_motion] = dec[~valid_motion]
        result["gaia_position_propagated_to_neowise_epoch_sep_arcsec"] = (
            angular_separation_arcsec(
                candidate["ra_deg"], candidate["dec_deg"],
                propagated_ra, propagated_dec,
            )
        )
        result["wise_median_epoch_year_approx"] = wise_epoch_year
        result = result.sort_values(
            "gaia_position_propagated_to_neowise_epoch_sep_arcsec"
            if "pmra" in result.columns else "gaia_epoch_separation_arcsec"
        ).reset_index(drop=True)
        all_rows.append(result)
        nearest = result.iloc[0]
        summary_rows.append({
            "candidate_id": candidate_id,
            "status": "COMPLETED",
            "query_mode": mode,
            "gaia_sources_within_15arcsec": int(len(result)),
            "nearest_source_id": str(nearest["source_id"]),
            "nearest_epoch2016_sep_arcsec": float(
                nearest["gaia_epoch_separation_arcsec"]
            ),
            "nearest_propagated_sep_arcsec": float(
                nearest["gaia_position_propagated_to_neowise_epoch_sep_arcsec"]
            ),
            "phot_variable_flag": (
                str(nearest["phot_variable_flag"])
                if "phot_variable_flag" in result.columns else None
            ),
            "best_class_name": (
                str(nearest["best_class_name"])
                if "best_class_name" in result.columns
                and pd.notna(nearest["best_class_name"]) else None
            ),
            "best_class_score": (
                float(nearest["best_class_score"])
                if "best_class_score" in result.columns
                and pd.notna(nearest["best_class_score"]) else None
            ),
            "nearest_g_mag": (
                float(nearest["phot_g_mean_mag"])
                if "phot_g_mean_mag" in result.columns
                and pd.notna(nearest["phot_g_mean_mag"]) else None
            ),
            "nearest_ruwe": (
                float(nearest["ruwe"])
                if "ruwe" in result.columns and pd.notna(nearest["ruwe"]) else None
            ),
            "wise_median_epoch_year_approx": wise_epoch_year,
        })

    combined = (
        pd.concat(all_rows, ignore_index=True)
        if all_rows else pd.DataFrame()
    )
    write_csv(combined, out_dir / "gaia_dr3_matches.csv")
    summary_frame = pd.DataFrame(summary_rows)
    write_csv(summary_frame, out_dir / "gaia_dr3_candidate_summary.csv")
    return {
        "status": (
            "COMPLETED" if any(row.get("status") == "COMPLETED" for row in summary_rows)
            else "PARTIAL_OR_UNAVAILABLE"
        ),
        "catalogue": "Gaia DR3",
        "service": GAIA_TAP,
        "search_radius_arcsec": 15.0,
        "candidate_summaries": summary_rows,
        "errors": errors,
        "interpretation": (
            "Gaia positions are at reference epoch 2016.0; proper-motion "
            "propagation is applied approximately to the median NEOWISE epoch. "
            "Gaia variability flags/classes are optical catalogue evidence, "
            "not a replacement for infrared time-series confirmation."
        ),
        "matches_file": "gaia_dr3_matches.csv",
        "summary_file": "gaia_dr3_candidate_summary.csv",
    }


def _ztf_summary_table(
    candidate_id: str,
    candidate: dict[str, float],
    frame: pd.DataFrame,
    collection: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if frame.empty:
        return pd.DataFrame(), {
            "candidate_id": candidate_id,
            "status": "NO_LIGHTCURVE_ROWS_IN_CONE",
            "collection": collection,
        }

    lower = {str(c).lower(): c for c in frame.columns}
    ra_col = lower.get("ra")
    dec_col = lower.get("dec")
    oid_col = lower.get("oid") or lower.get("objectid") or lower.get("object_id")
    mjd_col = lower.get("mjd")
    mag_col = lower.get("mag")
    magerr_col = lower.get("magerr") or lower.get("mag_err")
    band_col = lower.get("filtercode") or lower.get("fid") or lower.get("band")
    catflags_col = lower.get("catflags") or lower.get("catflags_mask")

    work = frame.copy()
    if oid_col is not None:
        work["_object_id"] = work[oid_col].astype(str)
    else:
        work["_object_id"] = "position_cone_result"

    if ra_col and dec_col:
        ra_values = pd.to_numeric(work[ra_col], errors="coerce").to_numpy(dtype=float)
        dec_values = pd.to_numeric(work[dec_col], errors="coerce").to_numpy(dtype=float)
        work["_separation_arcsec"] = angular_separation_arcsec(
            candidate["ra_deg"], candidate["dec_deg"], ra_values, dec_values
        )
    else:
        work["_separation_arcsec"] = np.nan

    if catflags_col:
        flags = pd.to_numeric(work[catflags_col], errors="coerce")
        work["_catflags"] = flags
        work["_strict_good"] = flags.fillna(-1).eq(0)
    else:
        work["_catflags"] = np.nan
        work["_strict_good"] = False

    if mjd_col:
        work["_mjd"] = pd.to_numeric(work[mjd_col], errors="coerce")
    else:
        work["_mjd"] = np.nan
    if mag_col:
        work["_mag"] = pd.to_numeric(work[mag_col], errors="coerce")
    else:
        work["_mag"] = np.nan
    if magerr_col:
        work["_magerr"] = pd.to_numeric(work[magerr_col], errors="coerce")
    else:
        work["_magerr"] = np.nan
    work["_band"] = work[band_col].astype(str) if band_col else "unknown"

    summary_rows = []
    for object_id, object_rows in work.groupby("_object_id", sort=False):
        usable = object_rows[
            np.isfinite(object_rows["_mjd"])
            & np.isfinite(object_rows["_mag"])
            & np.isfinite(object_rows["_magerr"])
            & (object_rows["_magerr"] > 0)
        ]
        good = usable.loc[usable["_strict_good"]]
        record = {
            "candidate_id": candidate_id,
            "ztf_object_id": object_id,
            "rows_returned": int(len(object_rows)),
            "photometry_rows": int(len(usable)),
            "strict_catflags_zero_rows": int(len(good)),
            "bands": ",".join(sorted(set(object_rows["_band"].dropna().astype(str)))),
            "nearest_separation_arcsec": (
                float(object_rows["_separation_arcsec"].min())
                if np.isfinite(object_rows["_separation_arcsec"]).any() else np.nan
            ),
            "mjd_min": float(usable["_mjd"].min()) if len(usable) else np.nan,
            "mjd_max": float(usable["_mjd"].max()) if len(usable) else np.nan,
            "median_mag_strict": float(good["_mag"].median()) if len(good) else np.nan,
            "mag_range_strict": (
                float(good["_mag"].max() - good["_mag"].min()) if len(good) else np.nan
            ),
            "status": "PHOTOMETRY_AVAILABLE" if len(usable) else "NO_VALID_PHOTOMETRY",
        }
        if len(good) >= 2:
            record["strict_mag_reduced_chi2"] = reduced_chi_square(
                good["_mag"].to_numpy(dtype=float),
                good["_magerr"].to_numpy(dtype=float),
            )
        else:
            record["strict_mag_reduced_chi2"] = np.nan
        summary_rows.append(record)
    return pd.DataFrame(summary_rows), {
        "candidate_id": candidate_id,
        "status": "COMPLETED",
        "collection": collection,
        "rows_returned": int(len(frame)),
        "objects_returned": int(work["_object_id"].nunique()),
        "summary": summary_rows,
    }


def run_ztf_check(out_dir: Path) -> dict[str, Any]:
    all_raw: list[pd.DataFrame] = []
    all_summaries: list[pd.DataFrame] = []
    candidate_results: list[dict[str, Any]] = []
    errors: list[str] = []

    for candidate_id, candidate in CANDIDATES.items():
        successful: pd.DataFrame | None = None
        used_collection = "default_public_collection"
        attempts = ("ztf_dr24", "ztf_dr23", "default")
        for collection in attempts:
            params = {
                "POS": (
                    f"CIRCLE {candidate['ra_deg']:.10f} "
                    f"{candidate['dec_deg']:.10f} {ZTF_RADIUS_DEG:.10f}"
                ),
                "BANDNAME": "g,r",
                "FORMAT": "CSV",
            }
            if collection != "default":
                params["COLLECTION"] = collection
            try:
                response = requests.get(
                    ZTF_LIGHTCURVE_API,
                    params=params,
                    timeout=(20, 150),
                    headers={"User-Agent": "NEOWISE-revalidation/2.0 research"},
                )
                response.raise_for_status()
                parsed = parse_delimited_response(response.text, sep=",")
                if len(parsed.columns) == 1 and len(parsed) == 0:
                    raise ValueError("ZTF returned an empty/unparseable table.")
                if not parsed.empty and not (
                    "mjd" in {str(c).lower() for c in parsed.columns}
                ):
                    raise ValueError(
                        "ZTF response did not contain an MJD column; "
                        f"columns={list(parsed.columns)}"
                    )
                successful = parsed
                used_collection = collection
                break
            except Exception as exc:
                errors.append(
                    f"ZTF {candidate_id} collection {collection}: "
                    f"{type(exc).__name__}: {exc}"
                )

        if successful is None:
            candidate_results.append({
                "candidate_id": candidate_id,
                "status": "UNAVAILABLE",
            })
            continue

        successful.insert(0, "candidate_id", candidate_id)
        all_raw.append(successful)
        summary, result = _ztf_summary_table(
            candidate_id, candidate, successful, used_collection
        )
        if len(summary):
            all_summaries.append(summary)
        candidate_results.append(result)

    raw = pd.concat(all_raw, ignore_index=True) if all_raw else pd.DataFrame()
    summary = (
        pd.concat(all_summaries, ignore_index=True)
        if all_summaries else pd.DataFrame()
    )
    write_csv(raw, out_dir / "ztf_lightcurve_raw.csv")
    write_csv(summary, out_dir / "ztf_lightcurve_summary.csv")
    return {
        "status": (
            "COMPLETED" if any(r.get("status") == "COMPLETED" for r in candidate_results)
            else "PARTIAL_OR_UNAVAILABLE"
        ),
        "catalogue": "Zwicky Transient Facility public light curves",
        "service": ZTF_LIGHTCURVE_API,
        "search_radius_arcsec": 5.0,
        "bands": ["g", "r"],
        "candidate_results": candidate_results,
        "errors": errors,
        "interpretation": (
            "ZTF is optical photometry, not mid-infrared confirmation. A positional "
            "non-detection or missing light curve does not disprove NEOWISE W1/W2 "
            "variability. Strict catflags==0 counts are reported separately."
        ),
        "raw_file": "ztf_lightcurve_raw.csv",
        "summary_file": "ztf_lightcurve_summary.csv",
    }


def _open_fits_bytes(content: bytes) -> tuple[np.ndarray, Any]:
    with fits.open(io.BytesIO(content), memmap=False) as hdus:
        for hdu in hdus:
            if hdu.data is not None:
                data = np.squeeze(np.asarray(hdu.data))
                header = hdu.header.copy()
                if data.ndim == 2:
                    return data, header
    raise ValueError("FITS file contained no two-dimensional image HDU.")


def _aperture_metrics(
    intensity: np.ndarray,
    uncertainty: np.ndarray,
    mask: np.ndarray,
    header: Any,
    ra_deg: float,
    dec_deg: float,
    catalog_mag: float,
) -> dict[str, Any]:
    wcs = WCS(header)
    x, y = wcs.world_to_pixel_values(float(ra_deg), float(dec_deg))
    yy, xx = np.indices(intensity.shape)
    radius = np.hypot(xx - float(x), yy - float(y))
    aperture_all = (radius <= APERTURE_RADIUS_PIX) & np.isfinite(intensity)
    annulus_all = (
        (radius >= ANNULUS_INNER_PIX)
        & (radius <= ANNULUS_OUTER_PIX)
        & np.isfinite(intensity)
    )
    if mask.shape != intensity.shape or uncertainty.shape != intensity.shape:
        raise ValueError("Intensity, uncertainty and mask cutouts have mismatched shapes.")

    mask_u32 = np.asarray(mask).astype(np.uint32, copy=False)
    recommended_bad = (
        np.bitwise_and(mask_u32, np.uint32(NEOWISE_RECOMMENDED_MASK_BITS)) != 0
    )
    clean_aperture = aperture_all & ~recommended_bad
    clean_annulus = annulus_all & ~recommended_bad
    if aperture_all.sum() < 5 or clean_aperture.sum() < 5 or clean_annulus.sum() < 20:
        raise ValueError("Too few unmasked pixels for the chosen aperture/annulus.")

    def measure(aperture: np.ndarray, annulus: np.ndarray) -> dict[str, float]:
        background_values = intensity[annulus]
        background = float(np.nanmedian(background_values))
        bg_sigma = float(
            1.4826 * np.nanmedian(np.abs(background_values - background))
        )
        values = intensity[aperture] - background
        flux = float(np.nansum(values))
        peak = float(np.nanmax(intensity[aperture]))
        unc = uncertainty[aperture]
        finite_unc = np.isfinite(unc) & (unc >= 0)
        if finite_unc.sum() == 0:
            flux_error = float("nan")
            snr = float("nan")
        else:
            variance = float(np.nansum(np.square(unc[finite_unc])))
            variance += int(aperture.sum()) * bg_sigma**2
            flux_error = math.sqrt(max(variance, 0.0))
            snr = flux / flux_error if flux_error > 0 else float("nan")
        zero_point = float(header.get("MAGZP", float("nan")))
        aperture_mag = (
            zero_point - 2.5 * math.log10(flux)
            if flux > 0 and math.isfinite(zero_point) else float("nan")
        )
        aperture_mag_error = (
            1.0857362047581296 * flux_error / flux
            if flux > 0 and math.isfinite(flux_error) else float("nan")
        )
        return {
            "background_dn": background,
            "background_robust_sigma_dn": bg_sigma,
            "aperture_flux_dn": flux,
            "aperture_flux_error_dn": flux_error,
            "aperture_snr": snr,
            "peak_dn": peak,
            "magzp": zero_point,
            "aperture_mag_uncorrected": aperture_mag,
            "aperture_mag_error_uncorrected": aperture_mag_error,
            "catalog_mag": float(catalog_mag),
            "aperture_minus_catalog_mag": (
                aperture_mag - float(catalog_mag)
                if math.isfinite(aperture_mag) and math.isfinite(catalog_mag)
                else float("nan")
            ),
        }

    raw = measure(aperture_all, annulus_all)
    clean = measure(clean_aperture, clean_annulus)
    return {
        "pixel_x": float(x),
        "pixel_y": float(y),
        "image_nx": int(intensity.shape[1]),
        "image_ny": int(intensity.shape[0]),
        "aperture_pixels_raw": int(aperture_all.sum()),
        "aperture_pixels_masked_recommended_bits": int(
            (aperture_all & recommended_bad).sum()
        ),
        "aperture_pixels_clean": int(clean_aperture.sum()),
        "annulus_pixels_clean": int(clean_annulus.sum()),
        "mask_recommended_fraction_in_aperture": float(
            (aperture_all & recommended_bad).sum() / max(aperture_all.sum(), 1)
        ),
        "raw": raw,
        "recommended_masked": clean,
        "header_magzpunc": header.get("MAGZPUNC"),
        "header_mjd_obs": header.get("MJD_OBS"),
        "header_moonsep_deg": header.get("MOONSEP"),
        "header_gcd2saa_deg": header.get("GCD2SAA"),
        "header_inevents": str(header.get("INEVENTS", "")),
    }


def run_image_remeasurement(
    candidate_detections: pd.DataFrame,
    out_dir: Path,
) -> dict[str, Any]:
    """Download and remeasure top four influential observations per candidate."""
    images_dir = out_dir / "neowise_fits_cutouts"
    images_dir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    manifest_rows: list[dict[str, Any]] = []
    measurement_rows: list[dict[str, Any]] = []
    fetch_errors: list[str] = []

    for candidate_id in CANDIDATES:
        group = candidate_detections.loc[
            candidate_detections["legacy_candidate_id"].astype(str) == candidate_id
        ].copy()
        if group.empty:
            fetch_errors.append(f"Candidate {candidate_id}: no observation rows.")
            continue
        for name in ("w1mpro", "w1sigmpro"):
            group[name] = pd.to_numeric(group[name], errors="coerce")
        group = group[
            np.isfinite(group["w1mpro"])
            & np.isfinite(group["w1sigmpro"])
            & (group["w1sigmpro"] > 0)
        ].copy()
        center = float(group["w1mpro"].median())
        group["legacy_influence_z"] = (
            (group["w1mpro"] - center) / group["w1sigmpro"]
        ).abs()
        group = group.sort_values(
            ["legacy_influence_z", "mjd"], ascending=[False, True]
        ).head(IMAGE_TOP_N_PER_CANDIDATE)

        for _, observation in group.iterrows():
            scan_id = str(observation["scan_id"]).strip()
            frame_number = int(observation["frame_num"])
            ra_deg = float(observation["ra"])
            dec_deg = float(observation["dec"])
            mjd = float(observation["mjd"])
            frame_id = f"{scan_id}{frame_number:03d}"
            scangrp = scan_id[-2:]
            for band in IMAGE_BANDS:
                catalog_mag = float(observation[f"w{band}mpro"])
                base_path = (
                    f"{scangrp}/{scan_id}/{frame_number:03d}/"
                    f"{frame_id}-w{band}-"
                )
                product_payload: dict[str, bytes] = {}
                product_urls: dict[str, str] = {}
                product_status: dict[str, str] = {}
                for product in IMAGE_PRODUCTS:
                    extension = ".fits" if product == "int" else ".fits.gz"
                    filename = f"{frame_id}-w{band}-{product}-1b{extension}"
                    url = (
                        f"{NEOWISE_IMAGE_ROOT}/{base_path}{product}-1b{extension}"
                    )
                    params = {
                        "center": f"{ra_deg:.8f},{dec_deg:.8f}",
                        "size": "101pix",
                    }
                    last_error: Exception | None = None
                    # The first campaign received intermittent HTTP 503 responses
                    # from IRSA IBE. Retry transient 5xx/429/network errors with
                    # bounded backoff; permanent 4xx and invalid payloads are not retried.
                    for attempt in range(4):
                        try:
                            response = session.get(
                                url,
                                params=params,
                                timeout=(20, 55),
                                headers={"User-Agent": "NEOWISE-revalidation/2.0 research"},
                            )
                            response.raise_for_status()
                            if not response.content.startswith((b"SIMPLE", b"\x1f\x8b")):
                                raise ValueError(
                                    "Response did not begin with FITS/gzip signature."
                                )
                            product_payload[product] = response.content
                            product_urls[product] = response.url
                            product_status[product] = "DOWNLOADED"
                            path = images_dir / filename
                            path.write_bytes(response.content)
                            manifest_rows.append({
                                "candidate_id": candidate_id,
                                "frame_id": frame_id,
                                "scan_id": scan_id,
                                "frame_num": frame_number,
                                "band": band,
                                "product": product,
                                "url": response.url,
                                "file": str(path.relative_to(out_dir)),
                                "size_bytes": len(response.content),
                                "sha256": sha256_bytes(response.content),
                                "status": "DOWNLOADED",
                                "attempts": attempt + 1,
                            })
                            last_error = None
                            break
                        except Exception as exc:
                            last_error = exc
                            status_code = getattr(
                                getattr(exc, "response", None), "status_code", None
                            )
                            retryable = (
                                status_code in {429, 500, 502, 503, 504}
                                or isinstance(
                                    exc,
                                    (requests.Timeout, requests.ConnectionError),
                                )
                            )
                            if retryable and attempt < 3:
                                time.sleep((1.5, 3.0, 6.0)[attempt])
                                continue
                            break
                    if product not in product_payload:
                        exc = last_error or RuntimeError("Unknown image download failure")
                        product_status[product] = f"ERROR:{type(exc).__name__}"
                        fetch_errors.append(
                            f"{candidate_id} {frame_id} W{band} {product}: "
                            f"{type(exc).__name__}: {exc}"
                        )
                        manifest_rows.append({
                            "candidate_id": candidate_id,
                            "frame_id": frame_id,
                            "scan_id": scan_id,
                            "frame_num": frame_number,
                            "band": band,
                            "product": product,
                            "url": url,
                            "status": f"ERROR:{type(exc).__name__}",
                            "error": str(exc)[:500],
                        })

                if len(product_payload) != 3:
                    measurement_rows.append({
                        "candidate_id": candidate_id,
                        "frame_id": frame_id,
                        "scan_id": scan_id,
                        "frame_num": frame_number,
                        "mjd": mjd,
                        "band": band,
                        "catalog_mag": catalog_mag,
                        "legacy_influence_z": float(observation["legacy_influence_z"]),
                        "status": "INCOMPLETE_IMAGE_SET",
                        **{f"{key}_status": value for key, value in product_status.items()},
                    })
                    continue
                try:
                    intensity, header = _open_fits_bytes(product_payload["int"])
                    uncertainty, _ = _open_fits_bytes(product_payload["unc"])
                    mask, _ = _open_fits_bytes(product_payload["msk"])
                    metrics = _aperture_metrics(
                        intensity, uncertainty, mask, header,
                        ra_deg, dec_deg, catalog_mag,
                    )
                    row = {
                        "candidate_id": candidate_id,
                        "frame_id": frame_id,
                        "scan_id": scan_id,
                        "frame_num": frame_number,
                        "mjd_catalog": mjd,
                        "band": band,
                        "catalog_mag": catalog_mag,
                        "legacy_influence_z": float(observation["legacy_influence_z"]),
                        "status": "REMEASURED",
                        "aperture_radius_pix": APERTURE_RADIUS_PIX,
                        "annulus_inner_pix": ANNULUS_INNER_PIX,
                        "annulus_outer_pix": ANNULUS_OUTER_PIX,
                        "aperture_correction_applied": False,
                        "mask_bits_excluded": NEOWISE_RECOMMENDED_MASK_BITS,
                        **{k: v for k, v in metrics.items() if k not in {"raw", "recommended_masked"}},
                    }
                    for prefix in ("raw", "recommended_masked"):
                        for key, value in metrics[prefix].items():
                            row[f"{prefix}_{key}"] = value
                    measurement_rows.append(row)
                except Exception as exc:
                    fetch_errors.append(
                        f"{candidate_id} {frame_id} W{band} photometry: "
                        f"{type(exc).__name__}: {exc}"
                    )
                    measurement_rows.append({
                        "candidate_id": candidate_id,
                        "frame_id": frame_id,
                        "scan_id": scan_id,
                        "frame_num": frame_number,
                        "mjd_catalog": mjd,
                        "band": band,
                        "catalog_mag": catalog_mag,
                        "status": f"PHOTOMETRY_ERROR:{type(exc).__name__}",
                        "error": str(exc)[:500],
                    })

    manifest = pd.DataFrame(manifest_rows)
    measurements = pd.DataFrame(measurement_rows)
    write_csv(manifest, out_dir / "neowise_image_download_manifest.csv")
    write_csv(measurements, out_dir / "neowise_image_remeasurements.csv")

    successful = (
        int(measurements["status"].eq("REMEASURED").sum())
        if "status" in measurements.columns else 0
    )
    downloaded = (
        int(manifest["status"].eq("DOWNLOADED").sum())
        if "status" in manifest.columns else 0
    )
    return {
        "status": "COMPLETED" if successful else "PARTIAL_OR_UNAVAILABLE",
        "archive": NEOWISE_IMAGE_ROOT,
        "source_rows_selected": int(IMAGE_TOP_N_PER_CANDIDATE * len(CANDIDATES)),
        "expected_image_sets": int(IMAGE_TOP_N_PER_CANDIDATE * len(CANDIDATES) * len(IMAGE_BANDS)),
        "image_products_per_set": list(IMAGE_PRODUCTS),
        "image_products_downloaded": downloaded,
        "remeasurements_completed": successful,
        "remeasurements_total_rows": int(len(measurements)),
        "measurement_file": "neowise_image_remeasurements.csv",
        "download_manifest": "neowise_image_download_manifest.csv",
        "download_directory": "neowise_fits_cutouts/",
        "errors": fetch_errors,
        "interpretation": (
            "Aperture magnitudes are uncorrected for aperture losses and are "
            "diagnostic comparisons with the archive profile-fit magnitudes, "
            "not replacements for the NEOWISE profile-fit catalogue values. "
            "Image masks are reported both before and after a conservative "
            "recommended bit selection."
        ),
    }


def run_injection_recovery(
    detections: pd.DataFrame,
    out_dir: Path,
) -> dict[str, Any]:
    """Inject controlled sinusoidal signals into cadence-preserving stable controls."""
    assigned, groups = build_groups(detections, GROUP_RADIUS_ARCSEC)
    stable = groups[
        (groups["n_observations"] >= MIN_OBSERVATIONS)
        & (groups["w1_reduced_chi2"] <= 1.2)
        & (groups["good_quality_fraction"] >= 0.90)
        & (groups["median_w1_snr"] >= 5.0)
    ].copy()
    group_data = {
        key: frame.copy()
        for key, frame in assigned.groupby("source_key", sort=False)
    }
    rng = np.random.default_rng(RANDOM_SEED)
    injection_rows: list[dict[str, Any]] = []
    null_rows: list[dict[str, Any]] = []
    stable_baselines: list[dict[str, Any]] = []

    for _, source in stable.iterrows():
        source_key = str(source["source_key"])
        lightcurve = group_data[source_key].copy()
        t = pd.to_numeric(lightcurve["mjd"], errors="coerce").to_numpy(dtype=float)
        x = pd.to_numeric(lightcurve["w1mpro"], errors="coerce").to_numpy(dtype=float)
        e = pd.to_numeric(lightcurve["w1sigmpro"], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(t) & np.isfinite(x) & np.isfinite(e) & (e > 0)
        t, x, e = t[valid], x[valid], e[valid]
        if len(x) < MIN_OBSERVATIONS:
            continue
        order = np.argsort(t)
        t, x, e = t[order], x[order], e[order]
        center = weighted_mean(x, e)
        base_chi2 = reduced_chi_square(x, e)
        residuals = x - center
        stable_baselines.append({
            "source_key": source_key,
            "n_observations": int(len(x)),
            "time_span_days": float(t.max() - t.min()),
            "baseline_w1_reduced_chi2": base_chi2,
            "good_quality_fraction": float(source["good_quality_fraction"]),
            "median_w1_snr": float(source["median_w1_snr"]),
        })

        # Residual bootstrap preserves each source's observation count and
        # sampling errors. It is conditional on this small internally selected
        # stable-control set, not an estimate for the full sky population.
        for replicate in range(NULL_BOOTSTRAP_REPLICATES):
            synthetic_residual = rng.choice(residuals, size=len(residuals), replace=True)
            simulated = center + synthetic_residual
            null_rows.append({
                "source_key": source_key,
                "replicate": replicate,
                "n_observations": int(len(x)),
                "reduced_chi2": reduced_chi_square(simulated, e),
                "flagged_at_chi2_ge_2": bool(
                    reduced_chi_square(simulated, e) >= VARIABILITY_THRESHOLD
                ),
            })

        phase_values = np.linspace(0.0, 2.0 * np.pi, INJECTION_PHASES, endpoint=False)
        t0 = float(t.min())
        for amplitude in INJECTION_AMPLITUDES_MAG:
            for period in INJECTION_PERIODS_DAYS:
                for phase_index, phase in enumerate(phase_values):
                    signal = amplitude * np.sin(
                        2.0 * np.pi * (t - t0) / period + phase
                    )
                    simulated = x + signal
                    chi2 = reduced_chi_square(simulated, e)
                    injection_rows.append({
                        "source_key": source_key,
                        "amplitude_semi_mag": amplitude,
                        "period_days": period,
                        "phase_index": phase_index,
                        "n_observations": int(len(x)),
                        "time_span_days": float(t.max() - t.min()),
                        "baseline_reduced_chi2": base_chi2,
                        "injected_reduced_chi2": chi2,
                        "recovered_at_chi2_ge_2": bool(
                            chi2 >= VARIABILITY_THRESHOLD
                        ),
                    })

    injection_frame = pd.DataFrame(injection_rows)
    null_frame = pd.DataFrame(null_rows)
    baseline_frame = pd.DataFrame(stable_baselines)
    write_csv(injection_frame, out_dir / "injection_recovery_trials.csv")
    write_csv(null_frame, out_dir / "residual_bootstrap_null_trials.csv")
    write_csv(baseline_frame, out_dir / "injection_stable_controls.csv")

    if len(injection_frame):
        recovery = (
            injection_frame.groupby(
                ["amplitude_semi_mag", "period_days"], as_index=False
            )
            .agg(
                trials=("recovered_at_chi2_ge_2", "size"),
                recovered=("recovered_at_chi2_ge_2", "sum"),
                recovery_fraction=("recovered_at_chi2_ge_2", "mean"),
                median_injected_chi2=("injected_reduced_chi2", "median"),
            )
        )
    else:
        recovery = pd.DataFrame()
    write_csv(recovery, out_dir / "injection_recovery_summary.csv")

    null_summary = (
        null_frame.groupby("source_key", as_index=False)
        .agg(
            bootstrap_replicates=("replicate", "size"),
            null_chi2_median=("reduced_chi2", "median"),
            null_chi2_95th_percentile=(
                "reduced_chi2", lambda s: float(np.nanpercentile(s, 95))
            ),
            conditional_flag_fraction=("flagged_at_chi2_ge_2", "mean"),
        )
        if len(null_frame) else pd.DataFrame()
    )
    write_csv(null_summary, out_dir / "residual_bootstrap_null_summary.csv")

    return {
        "status": "COMPLETED" if len(injection_frame) else "NO_STABLE_CONTROLS",
        "stable_controls_used": int(len(baseline_frame)),
        "injection_trials": int(len(injection_frame)),
        "null_bootstrap_trials": int(len(null_frame)),
        "amplitudes_semi_magnitude": list(INJECTION_AMPLITUDES_MAG),
        "periods_days": list(INJECTION_PERIODS_DAYS),
        "phases_per_amplitude_period": INJECTION_PHASES,
        "null_bootstrap_replicates_per_control": NULL_BOOTSTRAP_REPLICATES,
        "fixed_detection_threshold_reduced_chi2": VARIABILITY_THRESHOLD,
        "random_seed": RANDOM_SEED,
        "summary_file": "injection_recovery_summary.csv",
        "null_summary_file": "residual_bootstrap_null_summary.csv",
        "interpretation": (
            "Recovery fractions quantify sensitivity only for these sparse "
            "cadences, uncertainties, amplitudes and periods. Controls are "
            "selected internally by low EXP-005 chi-square and are not "
            "independently confirmed non-variable. Null bootstrap fractions "
            "are conditional diagnostics, not a global false-positive rate."
        ),
    }


def run_heldout_field(out_dir: Path) -> dict[str, Any]:
    """Apply frozen EXP-005 rules to a non-overlapping NEOWISE sky rectangle."""
    raw, selected_columns = retrieve_region(HELDOUT_BOUNDS)
    detections, counts = prepare_detections(raw, exclude_candidate_1=False)
    assigned, groups = build_groups(detections, GROUP_RADIUS_ARCSEC)
    groups["validator_positive_fixed_thresholds"] = (
        (groups["n_observations"] >= MIN_OBSERVATIONS)
        & (groups["median_w1_snr"] >= 5.0)
        & (groups["good_quality_fraction"] >= 0.90)
        & (groups["w1_reduced_chi2"] >= VARIABILITY_THRESHOLD)
    )
    groups["validator_stable_control_fixed_thresholds"] = (
        (groups["n_observations"] >= MIN_OBSERVATIONS)
        & (groups["median_w1_snr"] >= 5.0)
        & (groups["good_quality_fraction"] >= 0.90)
        & (groups["w1_reduced_chi2"] <= 1.2)
    )
    write_csv(groups, out_dir / "heldout_source_groups.csv")
    detections.to_parquet(out_dir / "heldout_analysis_detections.parquet", index=False)
    summary = {
        "status": "COMPLETED",
        "field_bounds_deg": {
            "ra_min": HELDOUT_BOUNDS[0],
            "ra_max": HELDOUT_BOUNDS[1],
            "dec_min": HELDOUT_BOUNDS[2],
            "dec_max": HELDOUT_BOUNDS[3],
        },
        "data_source": "same public NEOWISE year11 single-exposure source archive",
        "selected_columns": selected_columns,
        "population_counts": counts,
        "source_like_groups_ge_10": int(len(groups)),
        "fixed_threshold_variable_candidates": int(
            groups["validator_positive_fixed_thresholds"].sum()
        ),
        "fixed_threshold_stable_controls": int(
            groups["validator_stable_control_fixed_thresholds"].sum()
        ),
        "fixed_threshold_variable_fraction": float(
            groups["validator_positive_fixed_thresholds"].mean()
        ),
        "w1_reduced_chi2_median": float(groups["w1_reduced_chi2"].median()),
        "w1_reduced_chi2_95th_percentile": float(
            groups["w1_reduced_chi2"].quantile(0.95)
        ),
        "fixed_thresholds": {
            "min_observations": MIN_OBSERVATIONS,
            "min_median_w1_snr": 5.0,
            "min_good_quality_fraction": 0.90,
            "variable_reduced_chi2_min": VARIABILITY_THRESHOLD,
            "stable_reduced_chi2_max": 1.2,
            "grouping_radius_arcsec": GROUP_RADIUS_ARCSEC,
        },
        "files": [
            "heldout_source_groups.csv",
            "heldout_analysis_detections.parquet",
        ],
        "interpretation": (
            "This is spatial holdout validation using the same NEOWISE survey, "
            "archive, instrument and release. It tests transfer of fixed "
            "thresholds to a separate sky field; it is not an independent "
            "instrumental confirmation or a test set with external truth labels."
        ),
    }
    write_json(summary, out_dir / "heldout_field_summary.json")
    return summary


def _run_track_safely(
    track_name: str,
    function: Any,
    *args: Any,
) -> dict[str, Any]:
    """Keep a single external-service/analysis failure from suppressing other tracks."""
    try:
        return function(*args)
    except Exception as exc:
        return {
            "status": "FAILED",
            "track": track_name,
            "error_type": type(exc).__name__,
            "error": str(exc)[:2000],
            "interpretation": (
                "This track did not complete. Its result must not be inferred "
                "from empty output files or from the success of other tracks."
            ),
        }


def run(input_dir: Path, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_detections = pd.read_csv(input_dir / "candidate_detections.csv")
    baseline_groups = pd.read_csv(input_dir / "source_groups_eps1p5.csv")
    analysis_detections = pd.read_parquet(input_dir / "analysis_detections.parquet")

    # Each of the five tracks is isolated. A failure in one remote catalogue
    # or archive endpoint must not prevent the remaining scientific checks.
    gcvs = _run_track_safely(
        "GCVS positive control",
        run_gcvs_check,
        analysis_detections, baseline_groups, output_dir,
    )
    gaia = _run_track_safely(
        "Gaia DR3 cross-match", run_gaia_check,
        candidate_detections, output_dir,
    )
    ztf = _run_track_safely(
        "ZTF cross-match", run_ztf_check, output_dir,
    )
    images = _run_track_safely(
        "NEOWISE single-exposure image remeasurement",
        run_image_remeasurement,
        candidate_detections, output_dir,
    )
    injections = _run_track_safely(
        "Synthetic injection and null experiments",
        run_injection_recovery,
        analysis_detections, output_dir,
    )
    heldout = _run_track_safely(
        "Spatial held-out field", run_heldout_field, output_dir,
    )

    track_values = [gcvs, gaia, ztf, images, injections, heldout]
    failed_tracks = [
        item.get("track", item.get("status", "UNKNOWN"))
        for item in track_values if item.get("status") in {"FAILED", "UNAVAILABLE", "PARTIAL_OR_UNAVAILABLE"}
    ]
    summary = {
        "project": "Evidence-First Astronomical Anomaly Detection",
        "experiment": "VALIDATION-CAMPAIGN-1-TO-5",
        "status": (
            "COMPLETED_WITH_DOCUMENTED_LIMITATIONS"
            if not failed_tracks else "COMPLETED_WITH_TRACK_FAILURES"
        ),
        "failed_or_partial_tracks": failed_tracks,
        "campaign_tracks": {
            "1_gcvs_positive_control": gcvs,
            "2_gaia_and_ztf_external_crossmatches": {
                "gaia_dr3": gaia,
                "ztf": ztf,
            },
            "3_single_exposure_image_remeasurement": images,
            "4_synthetic_injection_and_null_tests": injections,
            "5_spatially_heldout_field": heldout,
        },
        "scientific_limits": [
            "A GCVS positional match is a catalogue cross-match, not a new discovery.",
            "Gaia and ZTF are optical; lack of optical variability does not disprove W1/W2 variability.",
            "Gaia DR3 positions are at epoch 2016.0; proper-motion propagation is approximate.",
            "Uncorrected aperture photometry is diagnostic and must not replace profile-fit photometry.",
            "The internal low-chi-square controls are not an externally labelled stable population.",
            "Injection recovery depends on the chosen amplitudes, periods, phase sampling, and cadence.",
            "The held-out sky field uses the same NEOWISE instrument and data release.",
            "No threshold was tuned using the held-out field.",
            "Failed or unavailable tracks are not interpreted as null scientific results.",
        ],
    }
    write_json(summary, output_dir / "validation_campaign_summary.json")

    manifest_rows = []
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "validation_campaign_manifest.json":
            data = path.read_bytes()
            manifest_rows.append({
                "path": str(path.relative_to(output_dir)),
                "size_bytes": len(data),
                "sha256": sha256_bytes(data),
            })
    write_json(
        {"experiment": "VALIDATION-CAMPAIGN-1-TO-5", "files": manifest_rows},
        output_dir / "validation_campaign_manifest.json",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run external cross-matches, image remeasurement, injection tests and spatial holdout."
    )
    parser.add_argument(
        "--input-dir", type=Path,
        default=Path("artifacts/neowise-revalidation"),
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("artifacts/validation-campaign"),
    )
    args = parser.parse_args()
    summary = run(args.input_dir, args.output_dir)
    for track, value in summary["campaign_tracks"].items():
        print(f"{track}: {value.get('status', 'SEE SUBTRACKS')}")
    print(f"Campaign output: {args.output_dir / 'validation_campaign_summary.json'}")


if __name__ == "__main__":
    main()
