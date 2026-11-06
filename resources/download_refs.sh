#!/bin/bash
# Download the reference files used by models/dnabert2 into $IMPRINT_BASE/{ref,variants}.
# Resumable (wget -c). On Bridges-2 login nodes, outside downloads run at only about 2 MB/s.
set -e
BASE=${IMPRINT_BASE:-/ocean/projects/cis260266p/rpandit1/imprintmethyl}
mkdir -p "$BASE/ref" "$BASE/variants"

# UCSC genomes and liftOver chain
cd "$BASE/ref"
wget -c -nv https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
wget -c -nv https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
wget -c -nv https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz

# Mouse Genomes Project REL-1505 per-strain SNPs (GRCm38/mm10; no mm9 release on the EBI server).
# 129X1/SvJ (the strain in the Xie et al. crosses) was not sequenced by MGP; 129S1/SvImJ is the closest proxy.
cd "$BASE/variants"
E=https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/strain_specific_vcfs
for s in CAST_EiJ 129S1_SvImJ; do
  for x in vcf.gz vcf.gz.tbi vcf.gz.md5; do
    wget -c -nv "$E/$s.mgp.v5.snps.dbSNP142.$x"
  done
  md5sum -c <(awk -v f="$s.mgp.v5.snps.dbSNP142.vcf.gz" '{print $1"  "f}' "$s.mgp.v5.snps.dbSNP142.vcf.gz.md5")
done

# Imprinted DMRs (Xie et al. 2012, Cell, Table S2), converted to CSV in this repo
cp "$(dirname "$0")/xie2012_imprinted_dmrs_mm9.csv" "$BASE/ref/" 2>/dev/null || true
