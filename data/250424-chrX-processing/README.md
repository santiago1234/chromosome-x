# Processing Chromosome X for Local Ancestry

This folder contains a simple Snakemake workflow to **prepare Chromosome X VCF files** for local ancestry analysis.

### Main steps:
- **Remove pseudoautosomal regions (PAR1 and PAR2)** based on GRCh38.
- **Make male samples diploid** (duplicate haploid genotypes).
- **Rename chromosome 23 to chrX** for consistency.

### Input:
- Raw VCF: `tgp_chrX.vcf.gz`
- PAR regions: `data/chrX_PAR_GRCh38.bed`
- Rename file: `data/rename.txt`
- Script to make males diploid: `scripts/convert_chromosome_X.py`

### Output:
- Final processed VCF: `chrX.vcf.gz`

### Quick pipeline overview:
1. **Exclude PAR regions** with `bcftools view`.
2. **Convert males to diploid** using a Python script.
3. **Rename chromosome** and compress the final VCF.

Intermediate files and logs are saved under `scratch/`.

### Requirements:
- `bcftools`
- `python`
