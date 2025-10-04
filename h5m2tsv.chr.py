#!/usr/bin/env python3

#
# h5m2csv.py: Convert HDF5 SNP matrix to CSV
#
# (c) 2021 by Joffrey Fitz (joffrey.fitz@tuebingen.mpg.de),
# Max Planck Institute for Developmental Biology,
# Tuebingen, Germany
#
# edited BY Frederique White

import sys
import h5py, numpy

input_file = sys.argv[1] 
chromosome = sys.argv[2] 
index = int(chromosome) - 1
label = f'chr{chromosome}'

f = h5py.File(input_file,'r')

# Array of tuples with start/stop indices for each chromosome
chr_regions = f['positions'].attrs['chr_regions']

start = chr_regions[index][0]
end = chr_regions[index][1]

# List of SNPs positions
positions = f['positions'][start:end]
# Matrix of SNPs genotypes
snp_chunk = f["snps"][start:end, :]

# Print header
print("#Chromosome", "Positon", "Count_zeros", "Count_ones"
	, "\t".join(f['accessions'][:].astype(str)), sep="\t" )

# Loop over all positions
#for pos in numpy.nditer(my_chr["positions"]):
with numpy.nditer(positions, flags=['multi_index']) as it:
	for pos in it:
		i = it.multi_index[0]

		# Get the corresponding SNPs for that position
		snps = snp_chunk[i]

		# Count 0s in snps
		cnt_zeros = numpy.count_nonzero(snps==0)

		# Count 1s in snp
		cnt_ones = numpy.count_nonzero(snps==1)

		print(label, pos, cnt_zeros, cnt_ones
			, "\t".join(snps.astype(str)), sep="\t")
	
