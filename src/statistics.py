"""Core uncertainty-aware statistics for astronomical source analysis."""

import numpy as np


def weighted_mean(values, errors):
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)

    valid = (
        np.isfinite(values)
        & np.isfinite(errors)
        & (errors > 0)
    )

    values = values[valid]
    errors = errors[valid]

    if len(values) == 0:
        return np.nan

    weights = 1.0 / errors**2
    return np.sum(weights * values) / np.sum(weights)


def reduced_chi_square(values, errors):
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)

    valid = (
        np.isfinite(values)
        & np.isfinite(errors)
        & (errors > 0)
    )

    values = values[valid]
    errors = errors[valid]

    n = len(values)
    if n < 2:
        return np.nan

    mean_value = weighted_mean(values, errors)
    chi2 = np.sum(((values - mean_value) / errors) ** 2)

    return chi2 / (n - 1)


def standardized_residuals(values, errors, mean=None):
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)

    if mean is None:
        mean = weighted_mean(values, errors)

    z = np.full(values.shape, np.nan, dtype=float)

    valid = (
        np.isfinite(values)
        & np.isfinite(errors)
        & (errors > 0)
    )

    z[valid] = (values[valid] - mean) / errors[valid]
    return z


def robust_mad_sigma(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return np.nan

    median = np.median(values)
    return 1.4826 * np.median(np.abs(values - median))
