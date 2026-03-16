
library(ggplot2)
library(qqman)

### args
args = commandArgs(trailingOnly = T)
gemma_prefix = args[1]
phenotype_file = args[2]

annotations_file = "annotations/annotations.len.tsv.gz"


### open GEMMA associations file
gemma_df = read.csv(paste0(gemma_prefix, ".assoc.txt"), sep="\t")

### remove dummy variables and missingness info
gemma_df$allele1 = gemma_df$allele0 = gemma_df$n_miss = NULL

### Bonferroni
alpha = 0.05
bonf_threshold = alpha / nrow(gemma_df)

### Benjamini-Hochberg
gemma_df$qval <- p.adjust(gemma_df$p_lrt, method = "BH")

cat(" * ", nrow(subset(gemma_df, p_lrt<0.05)), " SNPs at p<0.05\n")
cat(" * ", nrow(subset(gemma_df, p_lrt<0.01)), " SNPs at p<0.01\n")
cat(" * ", nrow(subset(gemma_df, p_lrt<0.001)), " SNPs at p<0.001\n")
cat(" * ", nrow(subset(gemma_df, qval<0.05)), " SNPs at q<0.05\n")

### Manhattan plots
#png(paste0(gemma_prefix, ".manhattan_Wald.png"), width = 1000, height = 400)
#manhattan(gemma_df, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline=-log10(bonf_threshold), genomewideline=F, main="Wald test")
#dev.off()

png(paste0(gemma_prefix,".manhattan_likelihoodratio.png"), width = 1000, height = 400)
manhattan(gemma_df, chr="chr", bp="ps", snp="rs", p="p_lrt", suggestiveline=-log10(bonf_threshold), genomewideline=F, main="likelihood ratio test")
dev.off()

#png(paste0(gemma_prefix, ".manhattan_score.png"), width = 1000, height = 400)
#manhattan(gemma_df, chr="chr", bp="ps", snp="rs", p="p_score", suggestiveline=-log10(bonf_threshold), genomewideline=F, main="Score test")
#dev.off()

### plot distributions
pdf(paste0(gemma_prefix, ".distribution.pdf"), width = 3, height = 3)
ggplot(gemma_df, aes(x=af)) + geom_histogram(bins=100) + xlab("MAF") + ylab("SNPs")
ggplot(gemma_df, aes(x=beta)) + geom_histogram(bins=100)+ xlab("beta") + ylab("SNPs")
ggplot(gemma_df, aes(x=logl_H1)) + geom_histogram(bins=100)+ xlab("log-likelihood H1") + ylab("SNPs")
ggplot(gemma_df, aes(x=l_remle)) + geom_histogram(bins=100)+ xlab("restricted maximum likelihood estimate of lambda") + ylab("SNPs")
ggplot(gemma_df, aes(x=l_mle)) + geom_histogram(bins=100)+ xlab("maximum likelihood estimate of the variance") + ylab("SNPs")
ggplot(gemma_df, aes(x=p_lrt)) + geom_histogram(bins=100)+ xlab("LRT p-value") + ylab("SNPs")
dev.off()

### Add SNPs annotations
signif_gemma_df = subset(gemma_df, qval<0.05)
write.table(signif_gemma_df$rs, file=paste0(gemma_prefix, "_snp_list.txt"), sep="\t", row.names=F, quote=F, col.names=F)
shell_command_extract_annotation = paste0("zgrep -wf ",gemma_prefix,"_snp_list.txt ",annotations_file," > " ,gemma_prefix,"_tmp.tsv")
system(shell_command_extract_annotation)
annotation_df = read.csv(paste0(gemma_prefix,"_tmp.tsv"), sep = "\t", header=F)
colnames(annotation_df) = c("rs","region",
	"L_gene","L_tss_distance","L_length","L_start","L_end","L_name","L_type","L_description","L_strand",
	"R_gene","R_tss_distance","R_length","R_start","R_end","R_name","R_type","R_description","R_strand",
	"subtype","ref","alt")


signif_gemma_df = merge(signif_gemma_df, annotation_df, by="rs", all.x=T)
signif_gemma_df = cbind(SNP=paste0(signif_gemma_df$rs,":",signif_gemma_df$ref,":",signif_gemma_df$alt), signif_gemma_df)
signif_gemma_df$rs = signif_gemma_df$ref = signif_gemma_df$alt = signif_gemma_df$chr = signif_gemma_df$ps = signif_gemma_df$p_wald = signif_gemma_df$p_score = NULL
write.table(signif_gemma_df, file=paste0(gemma_prefix, ".asso.FDR.tsv"), sep="\t", row.names=F, quote=F)





stop()

### List FDR significant SNPs
top_hits = subset(gemma_df, qval < 0.05)
top_hits = unique(top_hits[order(top_hits$qval), ]$rs)
write.table(top_hits, file=paste0(gemma_prefix, ".tophits.txt"), row.names=F, quote=F, col.names=F)

### use plink to extract genotypes of top hits
extract_geno_command = paste0("plink2 --bfile data/dataset_homo --extract ",gemma_prefix,".tophits.txt --recode A-transpose --out ", gemma_prefix,".tophits")
extract_geno_command2 = paste0("cut -f2,7- ", gemma_prefix, ".tophits.traw | sed 's/0_//g' > ", gemma_prefix, ".tophits.doses.tsv")
system(extract_geno_command)
system(extract_geno_command2)


### plot top hits with boxplots
pheno_df = read.csv(phenotype_file, sep = "\t", header=F)[, c(1,2,6)] # plink.fam
colnames(pheno_df) = c("fam","ID","estimate")
genotype_df = read.csv(paste0(gemma_prefix, ".tophits.doses.tsv"), sep = "\t", row.names = 1, check.names = F)
snpid = rownames(genotype_df)
genotype_df = data.frame(t(genotype_df))
colnames(genotype_df) = snpid

df = merge(pheno_df, genotype_df, by.x="ID", by.y="row.names")

top_top_hits = top_hits[1:5]
pdf(paste0(gemma_prefix, ".boxplot_top_hits.pdf"), width=3.5, height=5)
for (hit in top_top_hits){
	df[, hit] = df[, hit] + 2
	df[, hit] = ifelse(df[, hit]==4, 0, df[, hit])
	print(ggplot(df, aes(x = as.factor(df[, hit]), y = estimate, fill = as.factor(df[, hit]))) + 
		geom_boxplot(width=0.5, color="grey", outlier.shape = NA) + labs(x = hit, y = "Phenotype") + geom_jitter(width=0.3) + 
		theme(legend.position = "none", plot.title = element_text(color="black", size=10, face="bold")))
}
dev.off()

shell_command_rm_tmp_file = paste0("rm ", gemma_prefix,"_snp_list.txt ",gemma_prefix,"_tmp.tsv ", gemma_prefix,".tophits.log ", gemma_prefix,".tophits.traw ", gemma_prefix,".tophits.doses.tsv ")
system(shell_command_rm_tmp_file)
