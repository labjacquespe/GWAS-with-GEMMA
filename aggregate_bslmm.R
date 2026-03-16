
# Rscript aggregate_bslmm.R < prefix>
# by FW 

### args
gemma_prefix = commandArgs(trailingOnly = T)[1]
maf_file = "data/dataset_homo_stats.afreq"
phenotype_file = "data/dataset_homo.fam"
covariate_file = "../data/plant_size_standarized2.tsv"

hyp1_df = read.csv(paste0(gemma_prefix, "1", ".hyp.txt"), sep="\t")
hyp2_df = read.csv(paste0(gemma_prefix, "2", ".hyp.txt"), sep="\t")
hyp3_df = read.csv(paste0(gemma_prefix, "3", ".hyp.txt"), sep="\t")
aggregated_df = rbind(hyp1_df, hyp2_df, hyp3_df)
write.table(aggregated_df, file=paste0(gemma_prefix, "_aggregated.hyp.txt"), sep="\t", quote=F, row.names=F)

list_df = list(hyp1=hyp1_df, hyp2=hyp2_df, hyp3=hyp3_df, aggregated=aggregated_df)

get_summary = function(data){
	avg = round(mean(data), 3)
	median = round(median(data), 3)
	sd = round(sd(data),3)
	CI = round(quantile(data, probs=c(0.025, 0.975)), 3)
	return(list(avg=avg, median=median, sd=sd, CI=CI))
}

test_list = expand.grid(dataset=names(list_df), stat=c("pve","pge","n_gamma"))
results = as.data.frame(matrix(data=NA, nrow=nrow(test_list), ncol=7))
colnames(results) = c("stat","subset","mean","SD","median","CI2.5", "CI97.5")

for (i in 1:nrow(test_list)) {
	stat = as.character(test_list[i, "stat"])
	dataset_name = as.character(test_list[i, "dataset"])
	df = list_df[[dataset_name]]
	s = get_summary(df[, stat])
	results[i, ] = c(stat, dataset_name, s$avg, s$sd, s$median, s$CI[1], s$CI[2])
}

write.table(results, file=paste0(gemma_prefix, "_aggregated.summary.txt"), sep="\t", quote=F, row.names=F)
write.table(results, file="table1.tsv", sep="\t", quote=F, row.names=F)



param1_df = read.csv(paste0(gemma_prefix, "1", ".param.txt"), sep="\t")[, c(2,1,3,7)]
param2_df = read.csv(paste0(gemma_prefix, "2", ".param.txt"), sep="\t")[, c(2,7)]
param3_df = read.csv(paste0(gemma_prefix, "3", ".param.txt"), sep="\t")[, c(2,7)]

param1_df$n = param1_df$gamma * 100000
param2_df$n = param2_df$gamma * 100000
param3_df$n = param3_df$gamma * 100000

param_df = merge(param1_df, param2_df, by="rs", all=T)
param_df = merge(param_df, param3_df, by="rs", all=T)
colnames(param_df) = c("rs","chr","pos","pip1","gamma1","pip2","gamma2","pip3","gamma3")
param_df$PIP = rowSums(param_df[,c("gamma1","gamma2", "gamma3")])/300000
param_df$gamma1 = param_df$gamma2 = param_df$gamma3 = NULL

param_df = subset(param_df, PIP>0)

cat(">0 ", nrow(param_df), "\n")
cat(">0.0001 ", nrow(subset(param_df, PIP>0.0001)), "\n")
cat(">0.001 ", nrow(subset(param_df, PIP>0.001)), "\n")
cat(">0.01 ", nrow(subset(param_df, PIP>0.01)), "\n")

library(qqman)
#png(paste0(gemma_prefix,"_aggregated.manhattan_pip.png"), width = 900, height = 300)
pdf(paste0(gemma_prefix,"_aggregated.manhattan_pip.pdf"), width = 9, height = 3)
manhattan(param_df, chr="chr", bp="pos", snp="rs", p="PIP", suggestiveline=F, genomewideline=F, logp=F, ylab="Posterior inclusion probabilities (PIP)")
dev.off()

param_df = subset(param_df, PIP>0.001)

maf_df = read.csv(maf_file, sep="\t")[, c(2,6)]
colnames(maf_df) = c("rs", "MAF")
param_df = merge(param_df, maf_df, by="rs", all.x=T)
param_df = param_df[, c("chr", "pos", "rs","MAF","PIP")]

write.table(param_df, file=paste0(gemma_prefix,"_aggregated.param.txt"), sep="\t", row.names=F, quote=F)
write.table(param_df, file="table2.tsv", sep="\t", row.names=F, quote=F)




top_hits = subset(param_df, PIP >= 0.01)
top_hits = unique(top_hits[order(top_hits$PIP,decreasing=T), ]$rs)
write.table(top_hits, file=paste0(gemma_prefix, "_aggregated.tophits.txt"), row.names=F, quote=F, col.names=F)

### use plink to extract genotypes of top hits
extract_geno_command = paste0("plink2 --bfile data/dataset_homo --extract ",gemma_prefix,"_aggregated.tophits.txt --recode A-transpose --out ", gemma_prefix,"_aggregated.tophits")
extract_geno_command2 = paste0("cut -f2,7- ", gemma_prefix, "_aggregated.tophits.traw | sed 's/0_//g' > ", gemma_prefix, "_aggregated.tophits.doses.tsv")
system(extract_geno_command)
system(extract_geno_command2)

pheno_df = read.csv(phenotype_file, sep = "\t", header=F)[, c(2,6)] # plink.fam
colnames(pheno_df) = c("ID","phenotype")
cov_df = read.csv(covariate_file, sep = "\t")
colnames(cov_df) = c("ID","plant_size")
df = merge(pheno_df, cov_df, by="ID")

genotype_df = read.csv(paste0(gemma_prefix, "_aggregated.tophits.doses.tsv"), sep = "\t", row.names = 1, check.names = F)
snpid = rownames(genotype_df)
genotype_df = data.frame(t(genotype_df))
colnames(genotype_df) = snpid
genotype_df = genotype_df[, top_hits]

for (snp in top_hits){
	genotype_df[, snp] = genotype_df[, snp] + 2
	genotype_df[, snp] = ifelse(genotype_df[, snp]==4, 0, genotype_df[, snp])
	genotype_df[, snp] = ifelse(genotype_df[, snp]==2, 1, genotype_df[, snp])
}

df = merge(df, genotype_df, by.x="ID", by.y="row.names")
write.table(df, file=paste0(gemma_prefix,"_aggregated.genotypes.txt"), sep="\t", row.names=F, quote=F)
write.table(df, file="table3.tsv", sep="\t", row.names=F, quote=F)

shell_command_rm_tmp_file = paste0("rm ", gemma_prefix,"_aggregated.tophits.txt ", gemma_prefix,"_aggregated.tophits.log ", gemma_prefix,"_aggregated.tophits.traw ", gemma_prefix,"_aggregated.tophits.doses.tsv ")
system(shell_command_rm_tmp_file)

library(ggplot2)
pdf(paste0(gemma_prefix, "_aggregated.boxplot_top_hits.pdf"), width=3.5, height=5)
for (hit in top_hits){
	print(ggplot(df, aes(x = as.factor(df[, hit]), y = phenotype, fill = as.factor(df[, hit]))) + 
		geom_boxplot(width=0.5, color="grey", outlier.shape = NA) + labs(x = hit, y = "Phenotype") + geom_jitter(width=0.3) + 
		theme(legend.position = "none", plot.title = element_text(color="black", size=10, face="bold")))
}
dev.off()

shell_command_rm_tmp_file = paste0("rm ", gemma_prefix,"_aggregated.tophits.txt ", gemma_prefix,"_aggregated.tophits.log ", gemma_prefix,"_aggregated.tophits.traw ", gemma_prefix,"_aggregated.tophits.doses.tsv ")
system(shell_command_rm_tmp_file)

