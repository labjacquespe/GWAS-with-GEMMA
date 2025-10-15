#bslmm_summary.R

library(ggplot2)
#library(qqman)

### args
gemma_prefix = commandArgs(trailingOnly = T)[1]

#gamma = read.csv(paste0(gemma_prefix, ".gamma.txt"), sep="\t")

param = read.csv(paste0(gemma_prefix, ".param.txt"), sep="\t")
#param$index = rownames(param)
param$n_miss = NULL

### gamma column is (Posterior inclusion probabilities)
param = subset(param, gamma >= 0.001)

annotations = readRDS("annotations/SNP.EnsemblGenomes23.alleles.rds")
annotations$subtype = NULL
annotated_df = merge(param, annotations, by.x="rs", by.y="ID", all.x=T)
annotated_df = annotated_df[, c(1,7,8,6,4,5,9,10)]
annotated_df = annotated_df[order(annotated_df$gamma, decreasing=T),]
write.table(annotated_df, file=paste0(gemma_prefix, ".pip.txt"), sep="\t", row.names=F, quote=F)


#pdf(paste0(gemma_prefix, ".distribution.pdf"), width = 3, height = 3)
#ggplot(param, aes(x=gamma)) + geom_histogram(bins=100) + xlab("PIP") + ylab("SNPs")
#dev.off()


### Use Manhattan plot to visualize PIP
#lim = round(max(clean_pip_df$PIP)+0.1, 1)
#top = round(max(clean_pip_df$PIP)-0.2, 1)
#pdf(paste0(gemma_prefix,".manhattan_pip.pdf"), width = 7, height = 3)
#manhattan(clean_pip_df, chr="chr", bp="ps", snp="rs", p="PIP", suggestiveline=F, genomewideline=F, logp=F, ylab="PIP", ylim=c(0, lim), annotatePval=top)
#dev.off()



### Extract SNPs included in model with Proportion of variance explained (pve) > 0.9
hyp_df = read.csv(paste0(gemma_prefix, ".hyp.txt"), sep="\t")
hyp_df$iter = rownames(hyp_df)
cat("average Proportion of variance explained (PVE): ", mean(hyp_df$pve), "(SD ",sd(hyp_df$pve),")\n")
cat("average Proportion of genetic variance explained (PGE): ", mean(hyp_df$pge), "(SD ",sd(hyp_df$pge),")\n")
cat("average number of SNPs in models (n_gamma): ", mean(hyp_df$n_gamma), "(SD ",sd(hyp_df$n_gamma),")\n")

#pdf(paste0(gemma_prefix, ".diagnostic.pdf"), width = 3, height = 3)
png(paste0(gemma_prefix, ".diagnostic_pve.png"), width = 1000)
ggplot(hyp_df, aes(x = iter, y = pve)) + geom_point(alpha = 0.6) +
  labs(title = "PVE trace plot", x = "Iteration", y = "PVE") +
  theme_minimal()
dev.off()
 
png(paste0(gemma_prefix, ".diagnostic_pge.png"), width = 1000)
ggplot(hyp_df, aes(x = iter, y = pge)) + geom_point(alpha = 0.6) +
      labs(title = "PGE trace plot", x = "Iteration", y = "PGE") +
      theme_minimal()
dev.off()

png(paste0(gemma_prefix, ".diagnostic_n.png"), width = 1000)
ggplot(hyp_df, aes(x = iter, y = n_gamma)) +
  geom_point(alpha = 0.6) +
  theme_minimal()
dev.off()

#high_pve_iter = subset(hyp_df, pve > 0.9)

#snp_lists = list()
#for (i in 1:nrow(high_pve_iter)) {
#   iter = rownames(high_pve_iter)[i]
#   nsnp = high_pve_iter[i, "n_gamma"]
#   snp_lists[[iter]] = as.character(unlist(gamma[iter, 1:nsnp]))
#}
#param[snp_lists[[1]],]
