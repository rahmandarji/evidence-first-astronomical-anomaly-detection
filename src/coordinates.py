"""Coordinate utilities for astronomical source matching."""

import numpy as np


def parse_ra(value):
    try:
        s = str(value).strip()
        if not s or s.lower() in {"nan", "none"}:
            return np.nan
        import re
        s = re.sub(r"[hH]", ":", s)
        s = re.sub(r"[mM]", ":", s)
        s = re.sub(r"[sS]", "", s)
        s = re.sub(r"\s+", ":", s)
        s = re.sub(r":+", ":", s).strip(":")
        parts = s.split(":")
        if len(parts) == 3:
            h, m, sec = map(float, parts)
            return 15.0 * (h + m / 60.0 + sec / 3600.0)
        return float(s)
    except Exception:
        return np.nan


def parse_dec(value):
    try:
        s = str(value).strip()
        if not s or s.lower() in {"nan", "none"}:
            return np.nan
        import re
        sign = -1.0 if s.startswith("-") else 1.0
        s = s.lstrip("+-")
        s = re.sub(r"[dD]", ":", s)
        s = re.sub(r"[mM]", ":", s)
        s = re.sub(r"[sS]", "", s)
        s = re.sub(r"\s+", ":", s)
        s = re.sub(r":+", ":", s).strip(":")
        parts = s.split(":")
        if len(parts) == 3:
            d, m, sec = map(float, parts)
            return sign * (d + m / 60.0 + sec / 3600.0)
        return sign * abs(float(s))
    except Exception:
        return np.nan


def angular_distance_arcsec(ra1, dec1, ra2, dec2):
    ra1 = np.asarray(ra1, dtype=float)
    dec1 = np.asarray(dec1, dtype=float)
    ra2 = np.asarray(ra2, dtype=float)
    dec2 = np.asarray(dec2, dtype=float)
    mean_dec = np.deg2rad((dec1 + dec2) / 2.0)
    dra = (ra1 - ra2) * np.cos(mean_dec)
    ddec = dec1 - dec2
    return np.sqrt(dra**2 + ddec**2) * 3600.0
