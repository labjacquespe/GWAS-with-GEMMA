# GWAS in Arabidopsis thaliana with GEMMA

This github presents the steps to perform a GWAS using [GEMMA](https://github.com/genetics-statistics/GEMMA) (Genome-wide Efficient Mixed Model Association) on HPCs.

Here's the main steps included in this github:
1. Downloading genotypes and creating VCF file
2. Formatting genotypes
3. Calculating relatedness matrix
4. Running LLM
5. Running BSLMM
6. Annotating the results

Dependencies
- GEMMA
- Python 3
- PLINK2
- R
  
Python environment with the h5py package to read genotypes
(see *pyenv_requirements.txt*) 
```bash
ml python/3.11.5
virtualenv --no-download --clear pyenv
source pyenv/bin/activate
pip install --no-index h5py
```
PLINK and R
```bash
ml plink/2.00-20231024-avx2
ml gemma/0.98.5
ml r/4.4.0
ml mugqic/R_Bioconductor/4.3.2_3.18
```


### 1. Downloading genotypes from 1001 genomes (imputed matrix)
The imputed SNPs of 2029 _A. thalina_ lines can be downloaded at the [AraGWAS Catalog web site](https://aragwas.1001genomes.org/#/download-center).

1.1) Downloading (hdf5 format)
```bash
wget https://aragwas.1001genomes.org/api/genotypes/download
unzip download       # produce: GENOTYPES/4.hdf5
```

1.2) Extracting the SNPs from the 4.hdf5 file using a custom python script 
Extracting and merging the 5 chromosomes.
```bash
data_directory=data/GENOTYPES
sbatch extract_hdf5.sbatch 1:5 ${data_directory} # call h5m2tsv.chr.py

tail -n+2 ${data_directory}/genotypes.aragwas.chr2.tsv | cat ${data_directory}/genotypes.aragwas.chr1.tsv - > ${data_directory}/genotypes.aragwas.allchr.tsv
tail -n+2 ${data_directory}/genotypes.aragwas.chr3.tsv >> ${data_directory}/genotypes.aragwas.allchr.tsv
tail -n+2 ${data_directory}/genotypes.aragwas.chr4.tsv >> ${data_directory}/genotypes.aragwas.allchr.tsv
tail -n+2 ${data_directory}/genotypes.aragwas.chr5.tsv >> ${data_directory}/genotypes.aragwas.allchr.tsv
```

1.3) Converting to VCF format
The script takes a list of accessions to keep et keep SNPs with minor allele count > 0.
```bash
python create_vcf.py ${data_directory}/genotypes.aragwas.allchr.tsv data/included_samples.txt data/aragwas.genotypes_384.homozygote.vcf.gz
```


### 2) Formatting genotypes

2.1) Converting to PLINK 
```bash
plink2 --memory 10000 \
  --vcf data/aragwas.genotypes_384.homozygote.vcf.gz \
  --make-bed \
  --out data/aragwas.genotypes_384.homozygote
```
2.2) Adding phenotypes
```bash
plink2 --memory 10000 \
  --bfile data/aragwas.genotypes_384.homozygote \
  --make-bed \
  --pheno data/phenotypes.txt \
  --out data/dataset  
```

2.3) Looking at genotypes basic statistics
--freq will produce a report of minor allele frequencies
--hardy will produce a report of Hardy-Weinberg equilibrium
```bash
plink2 --bfile data/dataset --freq --hardy --out data/dataset_stats 
awk '$6>0.01' data/dataset_stats.afreq | wc -l
```


### 3) Calculating relatedness matrix
GEMMA provides relatedness matrix calculation. 
To estimate centered relationship, use the -gk1 parameters (details in section 4.4.1 of manual)

```bash
gemma -bfile data/dataset -gk 1 -r2 1.0 -o kinship -outdir output
```

### 4) Running LMM
Running GWAS with GEMMA, the frequentist test to use is optional (details in section 4.6 of manual):
* -lmm 1 performs Wald test
* -lmm 2 performs likelihood ratio test
* -lmm 3 performs score test
* -lmm 4 performs all the three tests

By default in GEMMA, SNPs with minor allele frequency below 1% will not be included in the analysis 
```bash
gemma -bfile data/dataset -k output/kinship.cXX.txt -lmm 4 -o lmm -outdir lmm/default
```

### 5) Running BSLMM
Fiting a Bayesian Sparse Linear Mixed Model (BSLMM) with GEMMA
Which model to fit (details in section 4.8 of manual):
* -bslmm 1 fits a standard linear BSLMM using MCMC
* -bslmm 2 fits a ridge regression/GBLUP with standard non-MCMC method
* -bslmm 3 fits a probit BSLMM using MCMC
```bash
gemma -bfile data/dataset -bslmm 1 -o bslmm -outdir bslmm/default
gemma -bfile data/dataset -bslmm 1 -w 2500000 -s 10000000 -rpace 100 -seed 1 -o chain1 -outdir bslmm/additional_iterations
```


### 6) Annotating the results

```bash
bash create_annotations.sh
```

Custom R script to produce Manhattan plots, boxplot of top hits and add SNP annotations

```bash
Rscript summary_lmm.R <output_prefix> <phenotypes>
Rscript summary_lmm.R lmm/default/lmm data/dataset.fam
Rscript summary_lmm.R lmm/adjusted_output/lmm data/dataset.fam
```
```bash
awk '{OFS="\t"; print $1,$3,$7,$14}' lmm/default/lmm.assoc.txt > table4.tsv
awk '{OFS="\t"; print $1,$3,$7,$14}' lmm/adjusted_output/lmm.assoc.txt > table5.tsv
```

Custom R script to compute posterior inclusion probabilities (PIP)
```bash
Rscript summary_bslmm.R bslmm/default
```
* You might need to run ```bash sed -i 's/\t$//g' bslmm/default.hyp.txt``` to remove tab at end of line

Custom R script to aggregate chains and compute posterior inclusion probabilities (PIP) across chains
```bash
Rscript aggregate_bslmm.R bslmm/additional_iterations/chain
```
