# GWAS in Arabidopsis thaliana with GEMMA

This github presents the steps to perform a GWAS using [GEMMA](https://github.com/genetics-statistics/GEMMA) (Genome-wide Efficient Mixed Model Association) on HPCs.

Here's the main steps included in this github:
1. Download of genotypes and create VCF file
2. Format the genotypes
3. Calculate relatedness matrix
4. Run LLM
5. Run BSLMM

Dependencies
- GEMMA
- Python 3
- PLINK2
- R 

Python environment with h5py package to read genotypes
```bash
ml python/3.11.5
virtualenv --no-download --clear gemma
source gemma/bin/activate
pip install --no-index h5py
```
PLINK and R
```bash
ml plink/2.00-20231024-avx2
ml gemma/0.98.5
ml r/4.4.0
ml mugqic/R_Bioconductor/4.3.2_3.18
```


### 1. Download the genotype from 1001 genomes (imputed matrix)
The imputed SNPs of 2029 _A. thalina_ lines can be downloaded at the [AraGWAS Catalog web site](https://aragwas.1001genomes.org/#/download-center).

1.1) Download file (hdf5 format)
```bash
wget https://aragwas.1001genomes.org/api/genotypes/download
unzip download       # produce: GENOTYPES/4.hdf5
```

1.2) Extract the SNPs from the GENOTYPES/4.hdf5 file using a custom python script 
Extract the 5 chromosomes then concat them.
```bash
python h5m2tsv.chr.py GENOTYPES/4.hdf5 1 > GENOTYPES/genotypes.aragwas.chr1.tsv
```

1.3) Transform to VCF format
The script takes a list of samples/lines to keep them remove SNPs with minor allele count of zero.
```bash
python create_vcf.py GENOTYPES/genotypes.aragwas.allchr.tsv data/included_samples.txt data/aragwas.genotypes_384.mac1.vcf.gz
```


### 2) Format the genotypes

2.1) Convert to PLINK 
```bash
plink2 --memory 10000 \
  --vcf data/aragwas.genotypes_384.mac1.vcf.gz \
  --make-bed \
  --out data/aragwas.genotypes_384.mac1
```
2.2) Add phenotype
```bash
plink2 --memory 10000 \
  --bfile data/aragwas.genotypes_384.mac1 \
  --make-bed \
  --pheno data/phenotypes.txt \
  --out data/dataset  
```

2.3) Look at genotypes basic statistics
--freq will produce a report of minor allele frequencies
--hardy will produce a report of Hardy-Weinberg equilibrium
```bash
plink2 --bfile data/dataset \
  --freq \
  --hardy \
  --out data/dataset_stats 
  
awk '$6>0.01' data/dataset_stats.afreq | wc -l
```


### 3) Calculate relatedness matrix
GEMMA provide the calculation of the relatedness matrix. 
To estimate centered relationship, use the -gk1 parameters (details in section 4.4.1 of manual)

```bash
gemma -bfile data/dataset -gk 1 -r2 1.0 -o kinship -outdir output
```


### 4) Run LLM
Running GWAS with GEMMA, the frequentist test to use is optional (details in section 4.6 of manual):
* -lmm 1 performs Wald test
* -lmm 2 performs likelihood ratio test
* -lmm 3 performs score test
* -lmm 4 performs all the three tests
```bash
gemma -bfile data/dataset -k output/kinship.cXX.txt -lmm 4 -o gemma_lmm4 -outdir output
```

Custom R script to produce Manhattan plots, boxplot of top hits and add SNP annotations
```bash
Rscript gwas_summary.R output/gemma_lmm4 data/phenotypes.txt
```


### 5) Run BSLMM
Fit a Bayesian Sparse Linear Mixed Model (BSLMM) with GEMMA
Which model to fit (details in section 4.8 of manual):
* -bslmm 1 fits a standard linear BSLMM using MCMC
* -bslmm 2 fits a ridge regression/GBLUP with standard non-MCMC method
* -bslmm 3 fits a probit BSLMM using MCMC
```bash
gemma -bfile data/dataset -bslmm 1 -o gemma_bslmm -outdir output
```

Custom R script to compute posterior inclusion probabilities (PIP)
```bash
Rscript bslmm_summary.R output/gemma_bslmm
```

