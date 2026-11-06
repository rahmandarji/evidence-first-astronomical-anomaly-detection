"""Rebuild the EXP-005 NEOWISE field and revalidate the two legacy candidates.

The original DBSCAN grouping and quality rule are preserved as a baseline.
Expanded IRSA quality filters are reported separately. Legacy group numbers
are reattached by sky coordinates, not by unstable DBSCAN label values.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow
from pyarrow import dataset as ds
from pyarrow.fs import S3FileSystem
from sklearn.cluster import DBSCAN

RA_CENTER, DEC_CENTER, RADIUS_DEG = 180.0, 0.0, 0.10
RA_MIN, RA_MAX = RA_CENTER - RADIUS_DEG, RA_CENTER + RADIUS_DEG
DEC_MIN, DEC_MAX = DEC_CENTER - RADIUS_DEG, DEC_CENTER + RADIUS_DEG
S3_BUCKET = "nasa-irsa-wise"
S3_PREFIX = "wise/neowiser/catalogs/p1bs_psd/healpix_k5/year11"
S3_METADATA = (
    f"{S3_BUCKET}/{S3_PREFIX}/"
    "neowiser-healpix_k5-year11.parquet/_metadata"
)
CANDIDATE_EXCLUSION = {
    "ra_deg": 179.936465, "dec_deg": -0.0897475, "radius_arcsec": 5.0
}
CANDIDATES = {
    "1978": {"ra_deg": 179.919108, "dec_deg": -0.053165},
    "1974": {"ra_deg": 179.928193, "dec_deg": -0.047839},
}
CANDIDATE_MATCH_MAX_ARCSEC = 5.0
GROUP_RADIUS_ARCSEC = 1.5
GROUPING_SENSITIVITY_ARCSEC = (1.0, 1.5, 2.0)
MIN_OBSERVATIONS = 10
CHI2_THRESHOLD = 2.0

REQUIRED_COLUMNS = [
    "ra", "dec", "mjd", "w1mpro", "w1sigmpro", "w1snr",
    "w2mpro", "w2sigmpro", "w2snr", "qual_frame", "ph_qual",
    "cc_flags", "scan_id", "frame_num",
]
OPTIONAL_COLUMNS = [
    "qi_fact", "saa_sep", "moon_masked", "w1rchi2", "w2rchi2",
    "source_id", "cntr",
]
IRSA_EXPANDED_QUALITY_FIELDS = ("qi_fact", "saa_sep", "moon_masked")


def _finite_pairs(values: Any, errors: Any) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(values, dtype=float)
    e = np.asarray(errors, dtype=float)
    if x.shape != e.shape:
        raise ValueError("values and errors must have identical shapes")
    valid = np.isfinite(x) & np.isfinite(e) & (e > 0)
    return x[valid], e[valid]


def weighted_mean(values: Any, errors: Any) -> float:
    """Inverse-variance weighted mean."""
    x, e = _finite_pairs(values, errors)
    if x.size == 0:
        return float("nan")
    weights = 1.0 / np.square(e)
    return float(np.sum(weights * x) / np.sum(weights))


def reduced_chi_square(values: Any, errors: Any) -> float:
    """Reduced chi-square around the inverse-variance weighted mean."""
    x, e = _finite_pairs(values, errors)
    if x.size < 2:
        return float("nan")
    center = weighted_mean(x, e)
    return float(np.sum(np.square((x - center) / e)) / (x.size - 1))


def arithmetic_center_reduced_chi_square(values: Any, errors: Any) -> float:
    """Legacy comparator: reduced chi-square around the arithmetic mean."""
    x, e = _finite_pairs(values, errors)
    if x.size < 2:
        return float("nan")
    center = float(np.mean(x))
    return float(np.sum(np.square((x - center) / e)) / (x.size - 1))


def angular_separation_arcsec(
    ra1_deg: float, dec1_deg: float, ra2_deg: Any, dec2_deg: Any
) -> np.ndarray:
    """Great-circle separation in arcseconds, stable across RA=0."""
    ra1, dec1 = math.radians(float(ra1_deg)), math.radians(float(dec1_deg))
    ra2 = np.radians(np.asarray(ra2_deg, dtype=float))
    dec2 = np.radians(np.asarray(dec2_deg, dtype=float))
    dra = (ra2 - ra1 + np.pi) % (2.0 * np.pi) - np.pi
    a = (
        np.sin((dec2 - dec1) / 2.0) ** 2
        + math.cos(dec1) * np.cos(dec2) * np.sin(dra / 2.0) ** 2
    )
    return 2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))) * (180.0 / np.pi) * 3600.0


def leave_one_out_chi_square(values: Any, errors: Any) -> np.ndarray:
    """Return weighted reduced chi-square after omitting each observation."""
    x, e = _finite_pairs(values, errors)
    if x.size < 3:
        return np.array([], dtype=float)
    output = []
    for omitted in range(x.size):
        keep = np.arange(x.size) != omitted
        output.append(reduced_chi_square(x[keep], e[keep]))
    return np.asarray(output, dtype=float)


def _package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "not-installed"


def retrieve_region() -> tuple[pd.DataFrame, list[str]]:
    """Retrieve the EXP-005 rectangle from the public IRSA S3 Parquet archive."""
    filesystem = S3FileSystem(region="us-west-2", anonymous=True)
    dataset = ds.parquet_dataset(
        S3_METADATA, filesystem=filesystem, partitioning="hive"
    )
    available = set(dataset.schema.names)
    missing = [name for name in REQUIRED_COLUMNS if name not in available]
    if missing:
        raise RuntimeError(
            "IRSA Parquet schema is missing required columns: " + ", ".join(missing)
        )
    columns = REQUIRED_COLUMNS + [
        name for name in OPTIONAL_COLUMNS
        if name in available and name not in REQUIRED_COLUMNS
    ]
    region_filter = (
        (ds.field("ra") >= RA_MIN)
        & (ds.field("ra") <= RA_MAX)
        & (ds.field("dec") >= DEC_MIN)
        & (ds.field("dec") <= DEC_MAX)
    )
    table = dataset.scanner(
        filter=region_filter, columns=columns, use_threads=True
    ).to_table()
    frame = table.to_pandas()
    if frame.empty:
        raise RuntimeError("IRSA returned zero rows for the documented field.")
    return frame, columns


def prepare_detections(raw: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Apply historical numeric, quality and Candidate #1 exclusion rules."""
    frame = raw.copy()
    numeric = [
        "ra", "dec", "mjd", "w1mpro", "w1sigmpro", "w1snr",
        "w2mpro", "w2sigmpro", "w2snr", "qual_frame", "frame_num",
        "qi_fact", "saa_sep", "w1rchi2", "w2rchi2",
    ]
    for name in numeric:
        if name in frame.columns:
            frame[name] = pd.to_numeric(frame[name], errors="coerce")
    for name in ("ph_qual", "cc_flags", "scan_id", "moon_masked"):
        if name in frame.columns:
            frame[name] = frame[name].fillna("").astype(str).str.strip()

    valid = (
        np.isfinite(frame["ra"]) & np.isfinite(frame["dec"])
        & np.isfinite(frame["mjd"]) & np.isfinite(frame["w1mpro"])
        & np.isfinite(frame["w1sigmpro"]) & (frame["w1sigmpro"] > 0)
    )
    usable = frame.loc[valid].copy().reset_index(drop=True)
    usable["original_good_quality"] = (
        (usable["qual_frame"] > 0)
        & (~usable["ph_qual"].str.contains("U", na=False))
        & (usable["cc_flags"] == "0000")
    )

    complete_expanded_fields = all(
        name in usable.columns for name in IRSA_EXPANDED_QUALITY_FIELDS
    )
    expanded = usable["original_good_quality"].copy()
    applied_rules: list[str] = []
    if "qi_fact" in usable.columns:
        expanded &= usable["qi_fact"].notna() & (usable["qi_fact"] > 0)
        applied_rules.append("qi_fact > 0")
    if "saa_sep" in usable.columns:
        expanded &= usable["saa_sep"].notna() & (usable["saa_sep"] > 0)
        applied_rules.append("saa_sep > 0")
    if "moon_masked" in usable.columns:
        expanded &= usable["moon_masked"].eq("00")
        applied_rules.append("moon_masked == '00'")
    usable["expanded_quality_available_filters"] = expanded.astype(bool)
    usable["expanded_quality_filter_complete"] = bool(complete_expanded_fields)

    exclusion_separation = angular_separation_arcsec(
        CANDIDATE_EXCLUSION["ra_deg"], CANDIDATE_EXCLUSION["dec_deg"],
        usable["ra"].to_numpy(), usable["dec"].to_numpy(),
    )
    exclusion_mask = exclusion_separation <= CANDIDATE_EXCLUSION["radius_arcsec"]
    excluded_count = int(exclusion_mask.sum())
    analysis = usable.loc[~exclusion_mask].copy().reset_index(drop=True)

    sort_columns = [
        name for name in ("ra", "dec", "mjd", "scan_id", "frame_num")
        if name in analysis.columns
    ]
    analysis = analysis.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)

    counts = {
        "raw_rows": int(len(raw)),
        "numerically_usable_rows": int(len(usable)),
        "original_good_quality_rows_before_exclusion": int(
            usable["original_good_quality"].sum()
        ),
        "expanded_quality_rows_before_exclusion": int(
            usable["expanded_quality_available_filters"].sum()
        ),
        "expanded_quality_fields_available": [
            name for name in IRSA_EXPANDED_QUALITY_FIELDS if name in usable.columns
        ],
        "expanded_quality_filter_complete": bool(complete_expanded_fields),
        "expanded_quality_rules_applied": applied_rules,
        "candidate_1_excluded_rows": excluded_count,
        "analysis_rows_after_candidate_1_exclusion": int(len(analysis)),
    }
    return analysis, counts


def _group_metrics(group: pd.DataFrame) -> dict[str, Any]:
    w1, w1err = _finite_pairs(
        group["w1mpro"].to_numpy(dtype=float),
        group["w1sigmpro"].to_numpy(dtype=float),
    )
    result: dict[str, Any] = {
        "n_observations": int(len(group)),
        "ra": float(group["ra"].median()),
        "dec": float(group["dec"].median()),
        "w1_weighted_mean": weighted_mean(w1, w1err),
        "w1_arithmetic_mean": float(np.mean(w1)),
        "w1_std": float(np.std(w1, ddof=1)) if len(w1) > 1 else float("nan"),
        "w1_median": float(np.median(w1)),
        "w1_range": float(np.ptp(w1)),
        "w1_reduced_chi2": reduced_chi_square(w1, w1err),
        "w1_arithmetic_center_reduced_chi2": arithmetic_center_reduced_chi_square(
            w1, w1err
        ),
        "w1_mad_sigma": float(
            1.4826 * np.median(np.abs(w1 - np.median(w1)))
        ),
        "median_w1_snr": float(
            np.nanmedian(pd.to_numeric(group["w1snr"], errors="coerce"))
        ),
        "median_w1_profile_fit_rchi2": (
            float(np.nanmedian(group["w1rchi2"]))
            if "w1rchi2" in group.columns else float("nan")
        ),
        "time_span_days": float(group["mjd"].max() - group["mjd"].min()),
        "good_quality_fraction": float(group["original_good_quality"].mean()),
        "expanded_quality_fraction": float(
            group["expanded_quality_available_filters"].mean()
        ),
        "expanded_quality_filter_complete": bool(
            group["expanded_quality_filter_complete"].all()
        ),
        "flagged_count": int((~group["original_good_quality"]).sum()),
        "n_scans": int(group["scan_id"].nunique()),
        "n_scan_frame_pairs": int(
            group[["scan_id", "frame_num"]].drop_duplicates().shape[0]
        ),
        "duplicate_scan_frame_pairs": int(
            group.groupby(["scan_id", "frame_num"], dropna=False).size().gt(1).sum()
        ),
    }

    w2_values = pd.to_numeric(group["w2mpro"], errors="coerce").to_numpy(dtype=float)
    w2_errors = pd.to_numeric(group["w2sigmpro"], errors="coerce").to_numpy(dtype=float)
    w2, w2err = _finite_pairs(w2_values, w2_errors)
    result["w2_n_observations"] = int(len(w2))
    result["w2_weighted_mean"] = weighted_mean(w2, w2err)
    result["w2_reduced_chi2"] = reduced_chi_square(w2, w2err)
    result["w2_median_snr"] = float(
        np.nanmedian(pd.to_numeric(group["w2snr"], errors="coerce"))
    )

    paired = group[
        np.isfinite(pd.to_numeric(group["w1mpro"], errors="coerce"))
        & np.isfinite(pd.to_numeric(group["w2mpro"], errors="coerce"))
    ]
    result["w1_w2_magnitude_correlation"] = (
        float(paired["w1mpro"].corr(paired["w2mpro"]))
        if len(paired) >= 3 else float("nan")
    )

    loo = leave_one_out_chi_square(w1, w1err)
    result["loo_n_tested"] = int(len(loo))
    result["loo_min_w1_reduced_chi2"] = (
        float(np.nanmin(loo)) if len(loo) else float("nan")
    )
    result["loo_max_w1_reduced_chi2"] = (
        float(np.nanmax(loo)) if len(loo) else float("nan")
    )
    result["loo_fraction_ge_2"] = (
        float(np.mean(loo >= CHI2_THRESHOLD)) if len(loo) else float("nan")
    )

    strict = group.loc[group["expanded_quality_available_filters"]]
    result["expanded_quality_n_observations"] = int(len(strict))
    result["expanded_quality_w1_reduced_chi2"] = reduced_chi_square(
        strict["w1mpro"].to_numpy(dtype=float),
        strict["w1sigmpro"].to_numpy(dtype=float),
    )
    result["expanded_quality_w1_weighted_mean"] = weighted_mean(
        strict["w1mpro"].to_numpy(dtype=float),
        strict["w1sigmpro"].to_numpy(dtype=float),
    )
    return result


def build_groups(
    detections: pd.DataFrame, eps_arcsec: float
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run original small-field DBSCAN, then assign deterministic local keys."""
    mean_dec_rad = np.deg2rad(float(detections["dec"].median()))
    coords = np.column_stack([
        (detections["ra"].to_numpy(dtype=float) - RA_CENTER)
        * np.cos(mean_dec_rad) * 3600.0,
        (detections["dec"].to_numpy(dtype=float) - DEC_CENTER) * 3600.0,
    ])
    labels = DBSCAN(
        eps=float(eps_arcsec), min_samples=1, metric="euclidean", n_jobs=1
    ).fit_predict(coords)

    assigned = detections.copy()
    assigned["_cluster_label"] = labels
    centroids = (
        assigned.groupby("_cluster_label", sort=False)
        .agg(ra=("ra", "median"), dec=("dec", "median"))
        .reset_index()
        .sort_values(["ra", "dec"], kind="mergesort")
        .reset_index(drop=True)
    )
    centroids["source_key"] = [
        f"eps{eps_arcsec:.1f}-S{index + 1:04d}"
        for index in range(len(centroids))
    ]
    key_map = dict(zip(centroids["_cluster_label"], centroids["source_key"]))
    assigned["source_key"] = assigned["_cluster_label"].map(key_map)

    rows = []
    for source_key, group in assigned.groupby("source_key", sort=False):
        if len(group) < MIN_OBSERVATIONS:
            continue
        row = {"source_key": source_key}
        row.update(_group_metrics(group))
        rows.append(row)
    groups = pd.DataFrame(rows)
    if groups.empty:
        raise RuntimeError(
            f"No groups with >= {MIN_OBSERVATIONS} detections at eps={eps_arcsec}."
        )
    groups = groups.sort_values(["ra", "dec"], kind="mergesort").reset_index(drop=True)
    return assigned.drop(columns=["_cluster_label"]), groups


def match_candidate(groups: pd.DataFrame, candidate: dict[str, float]) -> dict[str, Any]:
    """Match by saved coordinates, not legacy DBSCAN row-order IDs."""
    separations = angular_separation_arcsec(
        candidate["ra_deg"], candidate["dec_deg"],
        groups["ra"].to_numpy(dtype=float), groups["dec"].to_numpy(dtype=float),
    )
    index = int(np.argmin(separations))
    row = groups.iloc[index].to_dict()
    row["candidate_match_separation_arcsec"] = float(separations[index])
    row["candidate_match_status"] = (
        "MATCHED_WITHIN_5_ARCSEC"
        if separations[index] <= CANDIDATE_MATCH_MAX_ARCSEC
        else "NO_MATCH_WITHIN_5_ARCSEC"
    )
    return row


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.12g")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(output_dir: Path, write_raw: bool = True) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    raw, selected_columns = retrieve_region()
    detections, counts = prepare_detections(raw)

    if write_raw:
        raw.to_parquet(output_dir / "raw_region_detections.parquet", index=False)
    detections.to_parquet(output_dir / "analysis_detections.parquet", index=False)

    assigned, groups = build_groups(detections, GROUP_RADIUS_ARCSEC)
    write_csv(groups, output_dir / "source_groups_eps1p5.csv")

    candidate_rows = []
    observation_frames = []
    for legacy_id, coordinates in CANDIDATES.items():
        matched = match_candidate(groups, coordinates)
        matched["legacy_candidate_id"] = legacy_id
        candidate_rows.append(matched)
        if matched["candidate_match_status"] == "MATCHED_WITHIN_5_ARCSEC":
            selected = assigned.loc[assigned["source_key"] == matched["source_key"]].copy()
            selected.insert(0, "legacy_candidate_id", legacy_id)
            observation_frames.append(selected)

    write_csv(pd.DataFrame(candidate_rows), output_dir / "candidate_summary.csv")
    candidate_observations = (
        pd.concat(observation_frames, ignore_index=True)
        if observation_frames else pd.DataFrame()
    )
    write_csv(candidate_observations, output_dir / "candidate_detections.csv")

    sensitivity_rows = []
    for eps in GROUPING_SENSITIVITY_ARCSEC:
        _, eps_groups = build_groups(detections, eps)
        for legacy_id, coordinates in CANDIDATES.items():
            matched = match_candidate(eps_groups, coordinates)
            sensitivity_rows.append({
                "legacy_candidate_id": legacy_id,
                "eps_arcsec": eps,
                "n_source_like_groups": int(len(eps_groups)),
                "matched_source_key": matched["source_key"],
                "match_status": matched["candidate_match_status"],
                "match_separation_arcsec": matched["candidate_match_separation_arcsec"],
                "n_observations": matched["n_observations"],
                "w1_reduced_chi2": matched["w1_reduced_chi2"],
                "loo_min_w1_reduced_chi2": matched["loo_min_w1_reduced_chi2"],
                "expanded_quality_n_observations": matched["expanded_quality_n_observations"],
                "expanded_quality_w1_reduced_chi2": matched["expanded_quality_w1_reduced_chi2"],
            })
    write_csv(pd.DataFrame(sensitivity_rows), output_dir / "grouping_sensitivity.csv")

    summary = {
        "project": "Evidence-First Astronomical Anomaly Detection",
        "experiment": "REVALIDATION-V2",
        "status": "COMPUTATION_COMPLETE_PENDING_SCIENTIFIC_REVIEW",
        "data_source": {
            "archive": "NASA/IPAC Infrared Science Archive",
            "dataset": "NEOWISE Single-exposure Source Database",
            "release_partition": "year11",
            "s3_metadata_path": S3_METADATA,
            "region_deg": {
                "ra_min": RA_MIN, "ra_max": RA_MAX,
                "dec_min": DEC_MIN, "dec_max": DEC_MAX,
            },
            "selected_columns": selected_columns,
        },
        "population_counts": counts,
        "baseline_method": {
            "source_grouping": "scikit-learn DBSCAN",
            "eps_arcsec": GROUP_RADIUS_ARCSEC,
            "min_samples": 1,
            "minimum_observations": MIN_OBSERVATIONS,
            "quality_rule": "qual_frame > 0; no U in ph_qual; cc_flags == 0000",
            "variability_statistic": "weighted reduced chi-square around weighted mean",
        },
        "expanded_quality_method": {
            "additional_available_rules": counts["expanded_quality_rules_applied"],
            "all_three_irsa_fields_available": counts["expanded_quality_filter_complete"],
            "interpretation": "Additional sensitivity view, not a replacement baseline.",
        },
        "candidate_coordinate_matching": CANDIDATES,
        "grouping_sensitivity_arcsec": list(GROUPING_SENSITIVITY_ARCSEC),
        "candidates": candidate_rows,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "pyarrow": pyarrow.__version__,
            "scikit_learn": _package_version("scikit-learn"),
            "github_sha": os.environ.get("GITHUB_SHA", "local-unpinned"),
        },
        "scientific_limitations": [
            "Legacy group numbers are internal DBSCAN labels, not catalogue identifiers.",
            "Candidate IDs are reattached using stored median sky coordinates.",
            "DBSCAN can chain nearby detections and can merge or split apparent sources.",
            "IRSA quality filters reduce known risks but do not prove intrinsic variability.",
            "Reduced chi-square alone is not a calibrated probability or discovery claim.",
            "Independent confirmation requires external observations or catalogues.",
        ],
    }
    summary_path = output_dir / "run_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8"
    )

    artifacts = sorted(
        path for path in output_dir.iterdir()
        if path.is_file() and path.name != "manifest.json"
    )
    manifest = {
        "experiment": "REVALIDATION-V2",
        "files": [
            {
                "path": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in artifacts
        ],
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Retrieve and revalidate NEOWISE candidates 1978 and 1974."
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("artifacts/neowise-revalidation"),
    )
    parser.add_argument("--skip-raw", action="store_true")
    args = parser.parse_args()
    summary = run(args.output_dir, write_raw=not args.skip_raw)
    counts = summary["population_counts"]
    print("NEOWISE REVALIDATION V2")
    print(f"Raw detections: {counts['raw_rows']}")
    print(f"Usable W1 detections: {counts['numerically_usable_rows']}")
    print(f"Candidate #1 exclusions: {counts['candidate_1_excluded_rows']}")
    print(f"Analysis detections: {counts['analysis_rows_after_candidate_1_exclusion']}")
    print(f"Source-like groups: {len(pd.read_csv(args.output_dir / 'source_groups_eps1p5.csv'))}")
    print(pd.read_csv(args.output_dir / "candidate_summary.csv")[
        [
            "legacy_candidate_id", "candidate_match_status",
            "candidate_match_separation_arcsec", "n_observations",
            "w1_reduced_chi2", "loo_min_w1_reduced_chi2",
            "expanded_quality_n_observations", "expanded_quality_w1_reduced_chi2",
        ]
    ].to_string(index=False))
    print(f"Artifacts written to: {args.output_dir}")


if __name__ == "__main__":
    main()
