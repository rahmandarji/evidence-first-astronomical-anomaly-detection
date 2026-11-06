# NEOWISE candidate revalidation v2

This directory contains a separate, reproducible rebuild of the EXP-005 field.
It does not overwrite the historical notebooks, paper, or frozen results.

## Execution

The GitHub Actions workflow at .github/workflows/neowise-revalidation.yml runs
unit tests and then queries the public IRSA NEOWISE Parquet archive. It publishes
result files as a workflow artifact for 14 days.

For local execution with the pinned dependencies:

    python -m pip install -r revalidation/requirements.txt
    python revalidation/reproduce_neowise.py --output-dir artifacts/neowise-revalidation

No credentials or secrets are required. The data source is the public
nasa-irsa-wise S3 bucket, NEOWISE year11 Parquet metadata.

## Preserved baseline

The first pass intentionally retains the historical EXP-005 choices:

- RA 179.9 to 180.1 degrees; Dec -0.1 to +0.1 degrees.
- Valid W1 measurements require finite position/time/magnitude and positive
  W1 magnitude uncertainty.
- Original good-quality rule: qual_frame > 0, no U in ph_qual, and
  cc_flags == 0000.
- Candidate #1 is excluded within 5 arcseconds of its stored coordinates.
- DBSCAN with eps=1.5 arcsec, min_samples=1, and at least 10 detections per
  source-like group.
- Weighted reduced chi-square is calculated around the inverse-variance
  weighted mean.

## New diagnostics

1. Original EXP-005 quality rule.
2. Expanded quality view using any available IRSA fields among qi_fact > 0,
   saa_sep > 0, and moon_masked == 00.
3. Weighted-mean and arithmetic-mean chi-square results.
4. Leave-one-observation-out results for both candidate groups.
5. Candidate matching by stored sky coordinates, not unstable DBSCAN labels.
6. Grouping sensitivity at 1.0, 1.5, and 2.0 arcseconds.
7. Scan/frame coverage using the composite (scan_id, frame_num) key.

## Candidate identifiers

The values 1978 and 1974 are legacy internal group labels, not official
catalogue identifiers. This pipeline reattaches them by their saved median
coordinates:

- Legacy 1978: RA 179.919108 degrees, Dec -0.053165 degrees.
- Legacy 1974: RA 179.928193 degrees, Dec -0.047839 degrees.

A match is accepted only within 5 arcseconds. The resulting deterministic
source_key is local to this run, not a permanent astronomical source ID.

## Artifacts

- raw_region_detections.parquet: original field query before cleaning.
- analysis_detections.parquet: cleaned detections after Candidate #1 exclusion.
- source_groups_eps1p5.csv: source-group metrics for the baseline.
- candidate_summary.csv: candidate-level metrics and match status.
- candidate_detections.csv: observation rows for matched candidate groups.
- grouping_sensitivity.csv: candidate results under alternative grouping radii.
- run_summary.json: query, configuration, versions, counts, and limitations.
- manifest.json: SHA-256 checksums and sizes of the other result files.

## Interpretation limits

This is a reproducibility and revalidation computation, not a discovery claim.
An elevated reduced chi-square does not by itself establish intrinsic
astrophysical variability. Expanded quality filters reduce known observational
risks but cannot eliminate extraction, blending, cadence, calibration, or
instrumental explanations. Independent confirmation still requires external
observations or catalogues.
