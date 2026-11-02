#!/usr/bin/env python3
"""Summarize GEO's allele-resolved base-resolution methylC-seq calls."""

from __future__ import annotations

import argparse
import csv
import gzip
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path


SITE_HEADER = [
    "cross", "strand", "chrom", "pos_as_reported", "context",
    "129_methylated", "129_unmethylated", "129_coverage", "129_fraction",
    "Cast_methylated", "Cast_unmethylated", "Cast_coverage", "Cast_fraction",
    "maternal_methylated", "maternal_unmethylated", "maternal_coverage",
    "maternal_fraction", "paternal_methylated", "paternal_unmethylated",
    "paternal_coverage", "paternal_fraction",
    "unassigned_methylated", "unassigned_unmethylated",
]

SUMMARY_HEADER = [
    "cross", "strand", "context", "allele", "methylated_reads",
    "unmethylated_reads", "covered_site_observations", "methylation_fraction",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir", type=Path, default=Path("genetic_data/methylation_calls")
    )
    parser.add_argument("--output-dir", type=Path, default=Path("genetic_data/processed"))
    parser.add_argument(
        "--min-assigned-coverage",
        type=int,
        default=5,
        help="Keep CpG rows with at least this many 129+CAST reads.",
    )
    args = parser.parse_args()
    if args.min_assigned_coverage < 0:
        parser.error("--min-assigned-coverage must be non-negative")
    files = sorted(args.input_dir.glob("GSM7535*_mC_*.txt.gz"))
    if not files:
        raise FileNotFoundError(f"No processed methylation-call files found in {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    site_path = args.output_dir / "methylation_cpg_counts_mm9.tsv.gz"
    summary_path = args.output_dir / "methylation_context_summary.tsv"
    cpp_path = Path(__file__).with_name("summarize_methylation_calls.cpp")

    with tempfile.TemporaryDirectory(
        prefix="methylation-summary-", dir=args.output_dir
    ) as temporary_directory:
        temporary_path = Path(temporary_directory)
        executable = temporary_path / "summarize_methylation_calls"
        subprocess.run(
            [
                "g++",
                "-O3",
                "-std=c++17",
                str(cpp_path),
                "-o",
                str(executable),
            ],
            check=True,
        )
        processed_parts = []
        for path in files:
            name_parts = path.name.split("_")
            cross = name_parts[1]
            strand = name_parts[3].split(".")[0]
            site_part = temporary_path / f"{cross}_{strand}.sites.gz"
            summary_part = temporary_path / f"{cross}_{strand}.summary.tsv"
            command = (
                f"pigz -dc -- {shlex.quote(str(path))} | "
                f"{shlex.quote(str(executable))} "
                f"{shlex.quote(cross)} {shlex.quote(strand)} "
                f"{args.min_assigned_coverage} "
                f"{shlex.quote(str(summary_part))} | "
                f"pigz -1 > {shlex.quote(str(site_part))}"
            )
            subprocess.run(
                ["bash", "-o", "pipefail", "-c", command],
                check=True,
            )
            processed_parts.append((site_part, summary_part))

        temporary_site_path = temporary_path / site_path.name
        with gzip.open(
            temporary_site_path, "wt", newline="", compresslevel=1
        ) as site_output:
            csv.writer(site_output, delimiter="\t").writerow(SITE_HEADER)
        with temporary_site_path.open("ab") as site_output:
            for site_part, _ in processed_parts:
                with site_part.open("rb") as part:
                    shutil.copyfileobj(part, site_output)
        temporary_site_path.replace(site_path)

        with summary_path.open("w", newline="") as output:
            writer = csv.writer(output, delimiter="\t")
            writer.writerow(SUMMARY_HEADER)
            for _, summary_part in processed_parts:
                with summary_part.open(newline="") as summary_input:
                    writer.writerows(csv.reader(summary_input, delimiter="\t"))
    print(f"Wrote {site_path} and {summary_path}")


if __name__ == "__main__":
    main()
