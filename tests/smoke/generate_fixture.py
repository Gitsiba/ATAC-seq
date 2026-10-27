#!/usr/bin/env python3
"""Generate a deterministic, tiny paired-end ATAC-seq smoke-test fixture.

Synthetic reads are sampled from a narrow region of a random reference.
This fixture tests software wiring only; it is not biological benchmark data.
"""
from pathlib import Path
import csv
import random

ROOT = Path(__file__).resolve().parent
SEED = 27102026
REFERENCE_LENGTH = 20_000
READ_LENGTH = 75
FRAGMENT_LENGTH = 180
N_PAIRS = 300
ENRICHED_START = 8_000
ENRICHED_END = 9_000


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main():
    rng = random.Random(SEED)
    ROOT.mkdir(parents=True, exist_ok=True)
    reference = "".join(rng.choice("ACGT") for _ in range(REFERENCE_LENGTH))
    with (ROOT / "reference.fa").open("w") as handle:
        handle.write(">synthetic_chr1\n")
        for offset in range(0, len(reference), 80):
            handle.write(reference[offset:offset + 80] + "\n")

    r1_path = ROOT / "sample_A_R1.fastq"
    r2_path = ROOT / "sample_A_R2.fastq"
    max_start = ENRICHED_END - FRAGMENT_LENGTH
    with r1_path.open("w") as r1, r2_path.open("w") as r2:
        for index in range(N_PAIRS):
            start = rng.randint(ENRICHED_START, max_start)
            fragment = reference[start:start + FRAGMENT_LENGTH]
            read1 = fragment[:READ_LENGTH]
            read2 = reverse_complement(fragment[-READ_LENGTH:])
            name = f"@synthetic_{index:04d}"
            quality = "I" * READ_LENGTH
            r1.write(f"{name}/1\n{read1}\n+\n{quality}\n")
            r2.write(f"{name}/2\n{read2}\n+\n{quality}\n")

    with (ROOT / "samplesheet.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "fastq_1", "fastq_2"])
        writer.writerow(["sample_A", str(r1_path), str(r2_path)])

    print(f"Created test fixture in {ROOT}")


if __name__ == "__main__":
    main()
