# GEO methylation and imprinting data

The local GEO files are processed alignments, not raw sequencing reads. GEO's
`GSE30199_readme.txt` identifies the MethylC-seq coordinates as mm9. The RNA
tracks are 50-bp BED alignments from the reciprocal 129 × CAST frontal-cortex
crosses (`GSM828040` F1i and `GSM828041` F1r); each cross pooled three mice.
The RNA BED file list does not state its reference build, so this first pass
assumes mm9 based on the shared study/reference. Confirm the RNA build before
making definitive cross-assay overlaps.

## Reproducible first-pass processing

The mm9 RefSeq `refGene` annotation was downloaded from the UCSC mm9 database
to `genetic_data/reference/refGene.mm9.txt.gz`. To count the four allele-specific
RNA BED tracks over unique gene bodies:

```bash
python3 scripts/process_allelic_rna.py
```

Outputs are written to `genetic_data/processed/`:

- `rna_gene_body_counts_mm9.tsv`: gene coordinates, raw counts, and CPM using
  the number of records in each allele-specific BED as its library size. A BED
  read is counted only if it overlaps exactly one annotated gene body.
- `rna_bed_qc.tsv`: record counts and unique, ambiguous, and no-gene-body
  overlap totals for each input.
- `rna_reciprocal_allele_bias_mm9.tsv`: autosomal genes with at least 20
  uniquely assigned reads in each cross. The parent-origin and strain log2
  ratios are computed from library-normalized allele counts with a 0.5-count
  pseudocount. F1i is 129-maternal × CAST-paternal; F1r is CAST-maternal ×
  129-paternal.
- `rna_reciprocal_bias_scatter.png`: maternal-versus-paternal normalized
  allelic ratios in the two reciprocal crosses.

These are exploratory gene-body overlap counts, not exon-union gene counts,
transcript abundances, or replicate-level differential-expression tests.
Overlapping gene bodies are excluded from the count table's assigned reads;
sex chromosomes are retained in the count matrix but excluded from the
reciprocal parent-origin candidate table because the animals are male.

In this first pass, the four allele-specific BEDs contain 2.07–4.53 million
50-bp records each. About 77–78% of records overlap exactly one annotated gene
body; 10–14% overlap multiple gene bodies and are excluded from gene counts.
The candidate table contains 9,218 autosomal genes after the 20-read-per-cross
filter. Across these genes, reciprocal maternal-versus-paternal log2 ratios
have Pearson correlation -0.60, and 20% have the same direction in both
crosses. This suggests that strain-specific or technical allelic effects
dominate the genome-wide pattern; it is not evidence that most genes are
imprinted. A few known loci show the expected direction in this rough count
analysis (for example, Peg3/Peg10 paternal-like and Meg3 maternal-like), which
is a useful sanity check but not independent validation.

## Processed methylation calls

GEO also provides base-resolution calls as `*_mC_{crick,watson}.txt.gz`.
Unlike the BAM headers alone, these tables include explicit 129, CAST, and
unassigned methylated/unmethylated counts. The processing script is:

```bash
python3 scripts/summarize_methylation_calls.py
```

This streaming parser uses `g++` and `pigz` in addition to Python. It compiles
a small streaming parser, uses parallel gzip decompression, and writes the
multi-gigabyte output with low compression overhead.

It writes a read-weighted context summary and a CpG site-count table filtered
to sites with at least five assigned 129+CAST reads. The table retains
per-allele and maternal/paternal methylated and unmethylated counts, plus the
strand. Parent labels follow the reciprocal cross: F1i maternal=129 and
paternal=CAST; F1r maternal=CAST and paternal=129. The script processes every
matching Crick and Watson call file in the input directory and keeps the strand
label on every row.

Watson and Crick are the two complementary DNA strands. These files call
methylation on cytosines observed on each strand; at CpG sites, complementary
strand calls can describe the two cytosines in the same CpG dyad. The output
combines both files into one strand-labeled table, but does not add their counts
together or collapse them into one site. This avoids treating paired cytosines
as independent sites or merging coordinates incorrectly. For dyad-level
analyses, first verify GEO's coordinate convention and pair the complementary
calls using the appropriate offset.

`pos_as_reported` is intentionally not converted to BED coordinates; confirm
the source table's coordinate convention before joining it to BED tracks or
extracting sequence. If `chrL` occurs in the output, treat it as the lambda
spike-in, not a mouse chromosome; it can support conversion QC but must be
excluded from mouse-sequence modeling.

The four Crick/Watson files each contain about 525 million cytosine-call
records (about 2.10 billion total). With both strands included, the
five-assigned-read filter retained 6,928,034 F1i Crick, 6,968,834 F1i Watson,
6,893,853 F1r Crick, and 6,941,554 F1r Watson CpG observations. These
27,732,275 rows are strand- and cross-specific observations, not that many
unique genomic CpGs; complementary CpG cytosines have not been deduplicated.
Combining assigned 129+CAST reads across both strands, context methylation
fractions were F1i: 68.84% CG / 2.05% CHG / 2.77% CHH and F1r: 70.74% CG /
2.04% CHG / 2.80% CHH. The overall CG level still differs between pooled
crosses. These are descriptive QC summaries, not evidence of imprinting or
replicate-level significance.

The crosses are pooled samples, not independent per-mouse replicates. Treat
site-level and gene-level patterns as descriptive evidence, not independent
replicate-level significance tests or proof of causality. Standard bisulfite
sequencing also does not distinguish 5mC from 5hmC.

## Coordinate check and file note (added by Rutwik)

- **Coordinates are 1-based mm9.** Checked against the UCSC mm9 genome on chr1
  and chr19: watson `pos_as_reported` sits on the C of a CpG (19,867 of 20,000
  sampled sites, 99.3%), and crick sits on the G (26,554 of 26,886, 98.8%). So
  `crick_pos - 1` and `watson_pos` refer to the same CpG dyad.
- **Use `methylation_cpg_counts_mm9.tsv.gz`.** In the shared `processed.zip`,
  the uncompressed `methylation_cpg_counts_mm9.tsv` contains only the Crick rows
  (F1i 6,928,034, F1r 6,893,853). The `.tsv.gz` has all four Crick/Watson by
  cross tables (27,732,275 rows).
- **Imprinted DMRs** from Xie et al. 2012 Table S2 (55 DMRs, mm9; 41 maternally
  and 14 paternally methylated; 31 previously known) are in
  `resources/xie2012_imprinted_dmrs_mm9.csv`.
