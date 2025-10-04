
import sys
import gzip

infile = sys.argv[1]
included_samples_list = sys.argv[2]
outfile = sys.argv[3]

with open(included_samples_list) as f:
	keep_samples = set(line.strip() for line in f if line.strip())

with open(infile) as f_in, gzip.open(outfile, "wt") as f_out:
	header = f_in.readline().strip().split("\t")
	samples = header[4:]  # skip chromosome + position
	keep_indices = [i for i, s in enumerate(samples) if s in keep_samples]
	kept_samples = [samples[i] for i in keep_indices]

	# Write VCF header
	f_out.write("##fileformat=VCFv4.2\n")
	f_out.write("##source=HDF5_to_VCF\n")
	f_out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	f_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
	f_out.write("\t".join(kept_samples) + "\n")

	for line in f_in:
		parts = line.strip().split("\t")
		chrom, pos = parts[0].replace("chr", ""), parts[1]
		genotypes = [parts[4 + i] for i in keep_indices]

		gt_calls = []
		alt_count = 0
		for g in genotypes:
			if g == "0":
				gt_calls.append("0/0")
			elif g == "1":
				gt_calls.append("1/1")
				alt_count += 1
			else:  # missing or other
				gt_calls.append("./.")

			# Skip SNPs with no ALT alleles after filtering
		if alt_count == 0:
			continue
		
		# Dummy REF/ALT (N/*)
		ref, alt = "N", "*"
		snp_id = f"{chrom}:{pos}"

		f_out.write(f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gt_calls) + "\n")

