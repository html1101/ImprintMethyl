# ImprintMethyl

_Can genomic foundation models distinguish sequence-driven methylation from genomic imprinting?_

Question: Can a pretrained DNA model identify cis sequence differences associated with allele-specific methylation, and does it fail to predict methylation that depends on whether an allele came from the mother or father?

## Thing to keep in mind
- They use MethylC-Seq for measuring methylation data which can't distinguish between 5-methylcytosine (5mC) and 5-hydroxymethylcytosine (5hmC), 5hmC is prominent in the brain so important to note
- Methylation data uses mm9 reference genome (two different strains of it)
- Only was able to find 129 FASTA available for 129S1/SvlmJ, a proxy for the actual experiment's 129X1/SvJ strain.

## Data
1. `processed/ref/mm9-129.fa` - This contains the UCSC mm9 reference genome, modified specifically for the 129 strain.
2. `processed/ref/mm9-CAST.fa` - This contains the UCSC mm9 reference genome, modified specifically for the CAST strain.
3. 