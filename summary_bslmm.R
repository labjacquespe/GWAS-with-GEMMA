
# Rscript summary_bslmm.R < prefix>
# by FW 

library(ggplot2)
library(qqman)

### args
gemma_prefix = commandArgs(trailingOnly = T)[1]
annotations_file = "annotations/annotations.len.tsv.gz"

param = read.csv(paste0(gemma_prefix, ".param.txt"), sep="\t")
param$n_miss = NULL

### gamma column is (Posterior inclusion probabilities)
param = subset(param, gamma > 0)
cat("number of SNP with PIP > 0 ", nrow(param), "\n")

### Use Manhattan plot to visualize PIP
png(paste0(gemma_prefix,".manhattan_pip.png"), width = 1000, height = 400)
manhattan(param, chr="chr", bp="ps", snp="rs", p="gamma", suggestiveline=F, genomewideline=F, logp=F, ylab="PIP")
dev.off()


### extract SNP with PIP > 1%, annotate and write 
param = subset(param, gamma > 0.01)
write.table(param$rs, file=paste0(gemma_prefix, "_snp_list.txt"), sep="\t", row.names=F, quote=F, col.names=F)
tmp_file = paste0(gemma_prefix,"_tmp.tsv")
shell_command_extract_annotation = paste0("zgrep -wf ",gemma_prefix,"_snp_list.txt ",annotations_file," > " ,tmp_file)
system(shell_command_extract_annotation)
annotation_df = read.csv(tmp_file, sep = "\t", header=F)
colnames(annotation_df) = c("rs","region",
	"L_gene","L_tss_distance","L_length","L_start","L_end","L_name","L_type","L_description","L_strand",
	"R_gene","R_tss_distance","R_length","R_start","R_end","R_name","R_type","R_description","R_strand",
	"subtype","ref","alt")


param = merge(param, annotation_df, by="rs", all.x=T)
param = cbind(SNP=paste0(param$rs,":",param$ref,":",param$alt), param)
param$rs = param$ref = param$alt = param$chr = param$ps = NULL
write.table(param, file=paste0(gemma_prefix, ".PIP1pct.tsv"), sep="\t", row.names=F, quote=F)

shell_command_rm_tmp_file = paste0("rm ", gemma_prefix,"_snp_list.txt ",tmp_file)
system(shell_command_rm_tmp_file)

#param = subset(param, gamma > 0.01)
#cat("number of SNP with PIP > 1% ", nrow(param), "\n")
#annotations = readRDS("annotations/SNP.EnsemblGenomes23.alleles.rds")
#annotations$subtype = NULL
#annotated_df = merge(param, annotations, by.x="rs", by.y="ID", all.x=T)
#annotated_df = annotated_df[, c(1,7,8,6,4,5,9,10)]
#annotated_df = annotated_df[order(annotated_df$gamma, decreasing=T),]
#write.table(annotated_df, file=paste0(gemma_prefix, ".pip.txt"), sep="\t", row.names=F, quote=F)


### Extract SNPs included in model with Proportion of variance explained (pve) > 0.9
hyp_df = read.csv(paste0(gemma_prefix, ".hyp.txt"), sep="\t")
hyp_df$iter = rownames(hyp_df)
cat("average Proportion of variance explained (PVE): ", round(mean(hyp_df$pve),3), "(",round(sd(hyp_df$pve),3),") [",round(min(hyp_df$pve),3), "-", round(max(hyp_df$pve),3), "]\n")
cat("average Proportion of genetic variance explained (PGE): ", round(mean(hyp_df$pge),3), "(",round(sd(hyp_df$pge),3),") [",round(min(hyp_df$pge),3), "-", round(max(hyp_df$pge),3), "]\n")
cat("average number of SNPs in models (n_gamma): ", round(mean(hyp_df$n_gamma),3), "(",round(sd(hyp_df$n_gamma),3),") [",round(min(hyp_df$n_gamma),3), "-", round(max(hyp_df$n_gamma),3), "]\n")

### Trace plots
png(paste0(gemma_prefix, ".diagnostic_pve.png"), width = 2000)
ggplot(hyp_df, aes(x = iter, y = pve, group=1)) + geom_line()+ labs(title = "PVE trace plot", x = "Iteration", y = "PVE") + theme_minimal() + ylim(-0.1,1.1) 
dev.off()
 
png(paste0(gemma_prefix, ".diagnostic_pge.png"), width = 2000)
ggplot(hyp_df, aes(x = iter, y = pge, group=1)) + geom_line()+ labs(title = "PGE trace plot", x = "Iteration", y = "PGE") + theme_minimal() + ylim(-0.1,1.1) 
dev.off()

png(paste0(gemma_prefix, ".diagnostic_n.png"), width = 2000)
ggplot(hyp_df, aes(x = iter, y = n_gamma, group=1)) + geom_line()+ theme_minimal()
dev.off()
