#!/usr/bin/env python3
"""Summarize allele-specific RNA-seq BED alignments over mm9 gene bodies."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO


TRACKS = {
    "F1i_129": "GSM828040_F1i_RNA_Seq_129.bed",
    "F1i_CAST": "GSM828040_F1i_RNA_Seq_Cast.bed",
    "F1r_129": "GSM828041_F1r_RNA_Seq_129.bed",
    "F1r_CAST": "GSM828041_F1r_RNA_Seq_Cast.bed",
}
AUTOSOMES = {f"chr{number}" for number in range(1, 20)}


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def load_genes(path: Path):
    spans: dict[tuple[str, str], list[tuple[int, int, str]]] = defaultdict(list)
    for line_number, line in enumerate(open_text(path), start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 13:
            raise ValueError(f"{path}:{line_number}: expected at least 13 refGene fields")
        chrom, name = fields[2], fields[12]
        if not name or chrom in {"chrM", "chrL"}:
            continue
        try:
            start, end = int(fields[4]), int(fields[5])
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid transcript coordinates") from exc
        spans[(chrom, name)].append((start, end, fields[3]))

    intervals: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    gene_metadata = {}
    for (chrom, name), gene_spans in spans.items():
        gene_id = f"{chrom}:{name}"
        start = min(span[0] for span in gene_spans)
        end = max(span[1] for span in gene_spans)
        intervals[chrom].append((start, end, gene_id))
        strands = {span[2] for span in gene_spans}
        strand = strands.pop() if len(strands) == 1 else "."
        gene_metadata[gene_id] = (name, chrom, start, end, strand)
    for chrom in intervals:
        intervals[chrom].sort()
    return dict(intervals), gene_metadata


def make_index(
    intervals: dict[str, list[tuple[int, int, str]]],
) -> dict[str, tuple[list[int], list[int], list[str], list[int]]]:
    index = {}
    for chrom, entries in intervals.items():
        starts = [entry[0] for entry in entries]
        ends = [entry[1] for entry in entries]
        prefix_max_end = []
        current_max = 0
        for end in ends:
            current_max = max(current_max, end)
            prefix_max_end.append(current_max)
        index[chrom] = (starts, ends, [entry[2] for entry in entries], prefix_max_end)
    return index


def overlapping_genes(index, chrom: str, start: int, end: int) -> list[str]:
    if chrom not in index:
        return []
    starts, ends, genes, prefix_max_end = index[chrom]
    position = bisect.bisect_left(starts, end) - 1
    hits = []
    while position >= 0 and prefix_max_end[position] > start:
        if ends[position] > start:
            hits.append(genes[position])
        position -= 1
    return hits


def normalized_log2_ratio(
    numerator: int, numerator_library: int, denominator: int, denominator_library: int
) -> float:
    return math.log2(
        ((numerator + 0.5) / numerator_library)
        / ((denominator + 0.5) / denominator_library)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("genetic_data"))
    parser.add_argument(
        "--annotation",
        type=Path,
        default=Path("genetic_data/reference/refGene.mm9.txt.gz"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("genetic_data/processed"))
    parser.add_argument("--min-reads-per-cross", type=int, default=20)
    args = parser.parse_args()

    genes_by_chrom, gene_metadata = load_genes(args.annotation)
    index = make_index(genes_by_chrom)
    counts = {sample: Counter() for sample in TRACKS}
    qc = {}
    for sample, filename in TRACKS.items():
        path = args.data_dir / filename
        if not path.is_file():
            raise FileNotFoundError(f"Required RNA-seq BED is missing: {path}")
        metrics = Counter()
        chromosome_counts = Counter()
        with path.open("rt") as source:
            for line_number, line in enumerate(source, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    metrics["malformed"] += 1
                    continue
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    metrics["malformed"] += 1
                    continue
                if start < 0 or end <= start:
                    metrics["malformed"] += 1
                    continue
                metrics["reads"] += 1
                chromosome_counts[fields[0]] += 1
                hits = overlapping_genes(index, fields[0], start, end)
                if len(hits) == 1:
                    counts[sample][hits[0]] += 1
                    metrics["uniquely_assigned"] += 1
                elif len(hits) > 1:
                    metrics["ambiguous"] += 1
                else:
                    metrics["no_gene_body_overlap"] += 1
        metrics["chromosomes"] = len(chromosome_counts)
        metrics["malformed"] += 0
        qc[sample] = metrics
        print(
            f"{sample}: {metrics['reads']:,} reads; "
            f"{metrics['uniquely_assigned']:,} unique gene-body overlaps; "
            f"{metrics['ambiguous']:,} multi-gene overlaps"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    gene_ids = sorted(gene_metadata)
    library_sizes = {sample: qc[sample]["reads"] for sample in TRACKS}
    count_path = args.output_dir / "rna_gene_body_counts_mm9.tsv"
    with count_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(
            [
                "gene_id", "gene_symbol", "chrom", "start_mm9", "end_mm9", "strand",
                *TRACKS.keys(), *(f"{sample}_CPM" for sample in TRACKS),
            ]
        )
        for gene_id in gene_ids:
            symbol, chrom, start, end, strand = gene_metadata[gene_id]
            gene_counts = {sample: counts[sample][gene_id] for sample in TRACKS}
            writer.writerow(
                [
                    gene_id, symbol, chrom, start, end, strand,
                    *(gene_counts[sample] for sample in TRACKS),
                    *(
                        f"{gene_counts[sample] * 1_000_000 / library_sizes[sample]:.6g}"
                        for sample in TRACKS
                    ),
                ]
            )

    qc_path = args.output_dir / "rna_bed_qc.tsv"
    with qc_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(
            ["sample", "input_reads", "unique_gene_body_overlap", "multi_gene_overlap",
             "no_gene_body_overlap", "malformed_records", "chromosomes"]
        )
        for sample, metrics in qc.items():
            writer.writerow(
                [
                    sample,
                    metrics["reads"],
                    metrics["uniquely_assigned"],
                    metrics["ambiguous"],
                    metrics["no_gene_body_overlap"],
                    metrics["malformed"],
                    metrics["chromosomes"],
                ]
            )

    candidate_path = args.output_dir / "rna_reciprocal_allele_bias_mm9.tsv"
    with candidate_path.open("w", newline="") as output:
        columns = [
            "gene_id", "gene_symbol", "chrom", "start_mm9", "end_mm9",
            "F1i_129", "F1i_CAST", "F1r_129", "F1r_CAST",
            "F1i_129_fraction", "F1r_CAST_fraction",
            "F1i_maternal_log2_odds", "F1r_maternal_log2_odds",
            "mean_maternal_log2_odds", "mean_129_vs_CAST_log2_odds",
            "reciprocal_parent_direction_concordant",
        ]
        writer = csv.DictWriter(output, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for gene_id in gene_ids:
            symbol, chrom, start, end, _ = gene_metadata[gene_id]
            if chrom not in AUTOSOMES:
                continue
            f1i_129 = counts["F1i_129"][gene_id]
            f1i_cast = counts["F1i_CAST"][gene_id]
            f1r_129 = counts["F1r_129"][gene_id]
            f1r_cast = counts["F1r_CAST"][gene_id]
            if f1i_129 + f1i_cast < args.min_reads_per_cross:
                continue
            if f1r_129 + f1r_cast < args.min_reads_per_cross:
                continue
            i_ratio = normalized_log2_ratio(
                f1i_129,
                library_sizes["F1i_129"],
                f1i_cast,
                library_sizes["F1i_CAST"],
            )
            r_ratio = normalized_log2_ratio(
                f1r_129,
                library_sizes["F1r_129"],
                f1r_cast,
                library_sizes["F1r_CAST"],
            )
            f1i_parent = i_ratio
            f1r_parent = -r_ratio
            writer.writerow(
                {
                    "gene_id": gene_id,
                    "gene_symbol": symbol,
                    "chrom": chrom,
                    "start_mm9": start,
                    "end_mm9": end,
                    "F1i_129": f1i_129,
                    "F1i_CAST": f1i_cast,
                    "F1r_129": f1r_129,
                    "F1r_CAST": f1r_cast,
                    "F1i_129_fraction": f"{f1i_129 / (f1i_129 + f1i_cast):.6g}",
                    "F1r_CAST_fraction": f"{f1r_cast / (f1r_129 + f1r_cast):.6g}",
                    "F1i_maternal_log2_odds": f"{f1i_parent:.6g}",
                    "F1r_maternal_log2_odds": f"{f1r_parent:.6g}",
                    "mean_maternal_log2_odds": f"{(f1i_parent + f1r_parent) / 2:.6g}",
                    "mean_129_vs_CAST_log2_odds": f"{(i_ratio + r_ratio) / 2:.6g}",
                    "reciprocal_parent_direction_concordant": (
                        "yes" if f1i_parent * f1r_parent > 0 else "no"
                    ),
                }
            )
    print(f"Wrote {count_path}, {qc_path}, and {candidate_path}")


if __name__ == "__main__":
    main()
