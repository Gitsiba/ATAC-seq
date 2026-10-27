#!/usr/bin/env python3
"""Plot ATAC-seq fragment-length and peak-length distributions.

Real-data mode accepts a one-column fragment-length TSV and a MACS2
narrowPeak/BED-like file. Synthetic data are generated only with --demo.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_numeric_column(path: Path) -> np.ndarray:
    """Read one numeric value per row, tolerating one optional header."""
    values = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != 1:
                raise ValueError(
                    f"{path}:{line_number}: expected one numeric column, got {len(fields)}"
                )
            try:
                values.append(float(fields[0]))
            except ValueError:
                if line_number == 1 and not values:
                    continue
                raise ValueError(f"{path}:{line_number}: not a numeric value: {fields[0]!r}")
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        raise ValueError(f"No numeric observations found in {path}")
    if np.any(array <= 0):
        raise ValueError(f"Values in {path} must be positive")
    return array


def read_peak_lengths(path: Path) -> np.ndarray:
    """Read BED/narrowPeak columns 2–3 (0-based start, end)."""
    lengths = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least 3 BED columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: start/end must be integers") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval [{start}, {end})")
            lengths.append(end - start)
    array = np.asarray(lengths, dtype=int)
    if array.size == 0:
        raise ValueError(f"No valid intervals found in {path}")
    return array


def save_histogram(values: np.ndarray, title: str, xlabel: str, destination: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(values, bins=50, color="#5DADE2", edgecolor="white")
    ax.set(title=title, xlabel=xlabel, ylabel="Count")
    fig.tight_layout()
    fig.savefig(destination, dpi=160)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--insert-sizes", type=Path, help="One positive fragment length per row")
    parser.add_argument("--peaks", type=Path, help="BED/narrowPeak file; lengths use columns 2–3")
    parser.add_argument("--outdir", type=Path, default=Path("plots"), help="Output directory (default: plots)")
    parser.add_argument("--demo", action="store_true", help="Generate synthetic example plots, not experimental QC")
    args = parser.parse_args()

    if args.demo:
        if args.insert_sizes or args.peaks:
            parser.error("--demo cannot be combined with --insert-sizes or --peaks")
        rng = np.random.default_rng(42)
        insert_sizes = rng.lognormal(mean=4.5, sigma=0.35, size=3000).astype(int)
        peak_lengths = rng.lognormal(mean=5.4, sigma=0.55, size=2000).astype(int)
        label = " (synthetic demo data)"
    else:
        if not args.insert_sizes or not args.peaks:
            parser.error("provide both --insert-sizes and --peaks, or use --demo")
        insert_sizes = read_numeric_column(args.insert_sizes)
        peak_lengths = read_peak_lengths(args.peaks)
        label = ""

    args.outdir.mkdir(parents=True, exist_ok=True)
    save_histogram(insert_sizes, f"Insert size distribution{label}", "Fragment length (bp)",
                   args.outdir / "insert_size_distribution.png")
    save_histogram(peak_lengths, f"Peak length distribution{label}", "Peak length (bp)",
                   args.outdir / "peak_length_distribution.png")
    print(f"Saved QC plots to {args.outdir}")
    if args.demo:
        print("NOTE: plots use synthetic data and are not experimental QC results.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
