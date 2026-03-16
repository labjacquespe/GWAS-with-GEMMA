#!/bin/bash
#SBATCH --account=def-jacquesp
#SBATCH --mem-per-cpu=8G
#SBATCH --time=4:00:00
#SBATCH --chdir="./"

module load python/3.11.5
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r data/requirements.txt
pip install --no-index h5py

chrom=$1
workdir=$2
python h5m2tsv.chr.py ${workdir}/4.hdf5 ${chrom} > ${workdir}/genotypes.aragwas.chr${chrom}.tsv
