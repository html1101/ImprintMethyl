#!/bin/bash
# Download the reference files used by models/dnabert2 into $IMPRINT_BASE/{ref,variants}.
# Take the reference genome from mm10 and known SNPs for CAST and 129 strains, add them on top,
# and then map that from mm10 genome to mm9.
# use -c to make wget continue partially downloaded content

# On Bridges-2 login nodes, outside downloads run at only about 2 MB/s.
set -e
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE=${IMPRINT_BASE:-/ocean/projects/cis260266p/rpandit1/imprintmethyl}
mkdir -p "$BASE/ref" "$BASE/variants"

for tool in CrossMap bcftools samtools; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Required tool not found: $tool" >&2
    exit 1
  fi
done

# UCSC reference genomes and mm10-to-mm9 liftOver chain
cd "$BASE/ref"
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm9.over.chain.gz

# Mouse Genomes Project REL-1505 per-strain SNPs (GRCm38/mm10; no mm9 release on the EBI server).
# 129X1/SvJ (the strain in the Xie et al. crosses) was not sequenced by MGP; 129S1/SvImJ is the closest proxy.
cd "$BASE/variants"
E=https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs
for s in CAST_EiJ 129S1_SvImJ; do
  for x in vcf.gz vcf.gz.tbi vcf.gz.md5; do
    wget -c "$E/$s.mgp.v5.snps.dbSNP142.$x"
  done
  md5sum -c <(awk -v f="$s.mgp.v5.snps.dbSNP142.vcf.gz" '{print $1"  "f}' "$s.mgp.v5.snps.dbSNP142.vcf.gz.md5")
done

# The Mouse Genomes Project VCFs are in mm10 coordinates. Lift their SNPs back
# to mm9 before applying them to the mm9 reference used by the methylation data.
# 129S1/SvImJ is a proxy for the experiment's 129X1/SvJ.
REF_FA="$BASE/ref/mm9.fa"
if [[ ! -s "$REF_FA" ]]; then
  gzip -dc "$BASE/ref/mm9.fa.gz" > "$REF_FA"
fi
if [[ ! -s "$REF_FA.fai" || "$REF_FA" -nt "$REF_FA.fai" ]]; then
  samtools faidx "$REF_FA"
fi

TMP=$(mktemp -d "$BASE/ref/.strain-consensus.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

for s in CAST_EiJ 129S1_SvImJ; do
  VCF="$BASE/variants/$s.mgp.v5.snps.dbSNP142.vcf.gz"
  OUT="$BASE/ref/mm9_${s}.fa"
  CHAIN="$BASE/ref/mm10ToMm9.over.chain.gz"
  if [[ -s "$OUT" && ! "$VCF" -nt "$OUT" && ! "$REF_FA" -nt "$OUT" \
    && ! "$CHAIN" -nt "$OUT" ]]; then
    echo "Using existing $OUT"
    continue
  fi

  LIFTED_VCF="$TMP/$s.mm9.vcf"
  SORTED_VCF="$TMP/$s.mm9.sorted.vcf.gz"
  CrossMap vcf --chromid l "$CHAIN" "$VCF" "$REF_FA" "$LIFTED_VCF"
  bcftools view -f PASS -Oz -o "$TMP/$s.mm9.pass.vcf.gz" "$LIFTED_VCF"
  bcftools sort -Oz -o "$SORTED_VCF" "$TMP/$s.mm9.pass.vcf.gz"
  bcftools index "$SORTED_VCF"
  bcftools consensus --fasta-ref "$REF_FA" --sample "$s" --haplotype 1 \
    "$SORTED_VCF" > "$TMP/mm9_${s}.fa"
  mv "$TMP/mm9_${s}.fa" "$OUT"
done

# Imprinted DMRs (Xie et al. 2012, Cell, Table S2), converted to CSV in this repo
cp "$SCRIPT_DIR/xie2012_imprinted_dmrs_mm9.csv" "$BASE/ref/"
