# DNABERT-2 pipeline

Fine-tunes [DNABERT-2](https://github.com/MAGICS-LAB/DNABERT_2) (`zhihan1996/DNABERT-2-117M`) to predict per-allele CpG methylation from DNA sequence. It compares two conditions:

- **A: sequence only.** Input is a 1 kb window centred on the CpG.
- **B: sequence + parent of origin.** Same input plus a learned maternal/paternal embedding.

The gain of B over A is the predictive information carried by parent of origin beyond sequence. With pooled samples and one cross per direction, this is a predictive result, not a causal parental effect.

## Files

| File | What it does |
|---|---|
| `setup.sh` | Clones the `progen` conda env into `dnabert2`, installs `transformers==4.29.2`, `einops` and friends, and downloads the model weights into `$IMPRINT_BASE/hf_cache`. |
| `pilot.py`, `pilot.slurm` | Throughput test on random 1 kb sequences (no real data). |
| `finetune.py`, `finetune.slurm` | Builds the dataset from the processed methylation table, trains A and B, and writes `metrics.json`. |

All paths default to `/ocean/projects/cis260266p/rpandit1/imprintmethyl` and can be changed with `IMPRINT_BASE`. That folder is expected to contain:

- `processed/processed.zip`: Sarah's processed table. The script reads `methylation_cpg_counts_mm9.tsv.gz` from inside the zip. The uncompressed `.tsv` in the zip has only the crick strand.
- `ref/chromFa.tar.gz`: the UCSC mm9 genome.
- `ref/xie2012_imprinted_dmrs_mm9.csv`: a copy of `resources/xie2012_imprinted_dmrs_mm9.csv`.
- `dnabert2/`: these scripts.

`resources/download_refs.sh` downloads the UCSC mm9 and mm10 references, an
mm10-to-mm9 chain, and strain SNP VCFs. It lifts the SNPs from mm10 to mm9 and
builds `ref/mm9_CAST_EiJ.fa` and `ref/mm9_129S1_SvImJ.fa` against the mm9
reference, retaining PASS variants. The latter is a proxy for the experiment's
129X1/SvJ strain. The resulting FASTAs use mm9 coordinates and chromosome names,
so they align to the methylation table. `finetune.py` still uses the generic
mm9 reference for both alleles; changing it to select strain-specific sequence
is a separate step.

## Run on Bridges-2

```bash
bash setup.sh            # once
sbatch pilot.slurm       # optional speed check
sbatch finetune.slurm    # A vs B fine-tune; output in $IMPRINT_BASE/runs/run_<jobid>/metrics.json
```

Knobs, set as environment variables: `WIN` (window length, default 1000), `TRAIN_PCT` (percent of non-test CpGs sampled, default 5), `MAXTRAIN` (default 400000), `MAXTEST` (default 40000), `MINCOV` (minimum reads per allele, default 5), `EPOCHS` (default 1), `BATCH` (default 64).

## Data handling in `finetune.py`

- Uses autosomes only (chr1 to chr19), so chrX, chrY, chrM and the lambda spike-in are excluded.
- Merges strands per CpG. Coordinates are 1-based mm9: watson positions are the C of the CpG and crick positions are the G, so `crick_pos - 1` is the same CpG.
- Keeps CpG x cross pairs with at least `MINCOV` reads on both the maternal and the paternal allele. Each pair gives two examples (maternal, paternal).
- Splits:
  - **test:** all of chr19.
  - **dmr:** CpGs inside the 55 imprinted DMRs of Xie et al. 2012 (Table S2). These are held out of training.
  - **train:** a deterministic `TRAIN_PCT` sample of the remaining autosomes.
- Uses the mm9 reference sequence for both alleles. Strain-specific (129 / CAST) sequence is not used yet.

## Metrics

- MSE and Pearson r on chr19 and on the DMR set, against a predict-the-training-mean baseline.
- Correlation between predicted and observed maternal-minus-paternal differences. A predicts identical values for both alleles, so its difference is 0 by construction.
- At the DMR level: how many DMRs show the published methylated parent in the data, and how many in condition B's predictions.

## Pilot throughput

Measured on Bridges-2 with bf16, batch 32 and 1 kb windows (222 tokens, about 4.5 bases per token):

| GPU | Step time | Examples/s | Peak memory |
|---|---|---|---|
| H100 80 GB | 50 ms | about 640 | 7.1 GB |
| L40S 48 GB | 116 ms | about 277 | 7.0 GB |

DNABERT-2's Triton flash-attention kernel fails to compile on this stack (torch 2.7, CUDA 12.8). The scripts set `sys.modules["triton"] = None` so the model falls back to standard PyTorch attention. At about 640 examples/s this sustains roughly 1e14 FLOP/s on an H100, about 10% of peak.

## Results

First-pass fine-tune: pending (job running on Bridges-2).

## Known limitations

- Both alleles use the mm9 reference sequence, so sequence-dependent allelic effects are invisible to the model.
- The samples are pooled frontal cortex, one cross per direction, so there are no biological replicates.
- Bisulfite sequencing does not separate 5mC from 5hmC.
- Only one held-out test chromosome is used for now.
