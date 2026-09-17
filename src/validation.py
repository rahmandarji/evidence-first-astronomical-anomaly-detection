"""Uncertainty-aware variability validation utilities."""

import numpy as np


def calculate_validator_metrics(
    magnitudes,
    errors,
    quality_good=None,
):
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)

    valid = (
        np.isfinite(magnitudes)
        & np.isfinite(errors)
        & (errors > 0)
    )

    magnitudes = magnitudes[valid]
    errors = errors[valid]

    n = len(magnitudes)

    if n == 0:
        return {
            "n_obs": 0,
            "median_snr_proxy": np.nan,
            "median_w1_err": np.nan,
            "w1_mean": np.nan,
            "w1_std": np.nan,
            "w1_range": np.nan,
            "reduced_chi2": np.nan,
            "quality_fraction": np.nan,
        }

    mean_mag = np.mean(magnitudes)
    chi2 = np.sum(((magnitudes - mean_mag) / errors) ** 2)
    dof = max(n - 1, 1)

    if n > 1:
        w1_std = np.std(magnitudes, ddof=1)
    else:
        w1_std = 0.0

    if quality_good is not None:
        quality_good = np.asarray(quality_good, dtype=float)
        quality_fraction = float(np.mean(quality_good))
    else:
        quality_fraction = np.nan

    return {
        "n_obs": int(n),
        "median_snr_proxy": float(np.median(1.0 / errors)),
        "median_w1_err": float(np.median(errors)),
        "w1_mean": float(mean_mag),
        "w1_std": float(w1_std),
        "w1_range": float(np.max(magnitudes) - np.min(magnitudes)),
        "reduced_chi2": float(chi2 / dof),
        "quality_fraction": quality_fraction,
    }


def classify_variability(
    metrics,
    min_observations=10,
    min_snr_proxy=5.0,
    reduced_chi2_threshold=2.0,
    good_quality_fraction=0.90,
):
    variable_signal = (
        np.isfinite(metrics["reduced_chi2"])
        and metrics["reduced_chi2"] >= reduced_chi2_threshold
    )

    enough_observations = metrics["n_obs"] >= min_observations

    sufficient_snr = (
        np.isfinite(metrics["median_snr_proxy"])
        and metrics["median_snr_proxy"] >= min_snr_proxy
    )

    quality = metrics["quality_fraction"]

    good_quality = (
        np.isnan(quality)
        or quality >= good_quality_fraction
    )

    if variable_signal and enough_observations and sufficient_snr:
        label = (
            "VARIABLE_HIGH_QUALITY"
            if good_quality
            else "VARIABLE"
        )
    elif enough_observations and sufficient_snr:
        label = "NOT_VARIABLE_BY_CURRENT_GATE"
    else:
        label = "INSUFFICIENT_DATA"

    return {
        "validator_class": label,
        "variable_signal": bool(variable_signal),
        "enough_observations": bool(enough_observations),
        "sufficient_snr": bool(sufficient_snr),
        "good_quality": bool(good_quality),
    }
