#bslmm_summary.R

library(ggplot2)
library(qqman)

### args
gemma_prefix = commandArgs(trailingOnly = T)[1]

gamma = read.csv(paste0(gemma_prefix, ".gamma.txt"), sep="\t")

param = read.csv(paste0(gemma_prefix, ".param.txt"), sep="\t")
param$index = rownames(param)
param$n_miss = NULL

### compute PIP (Posterior inclusion probabilities)
vals = unlist(gamma)
vals = vals[vals != 0]
tab = table(vals)
n_iter = nrow(gamma)

pip_df = data.frame(snp_index = as.integer(names(tab)), count = as.integer(tab), PIP = as.integer(tab) / n_iter)
clean_pip_df = merge(param, pip_df, by.x="index", by.y="snp_index")
clean_pip_df$index=NULL
write.table(clean_pip_df, file=paste0(gemma_prefix, ".pip.txt"), sep="\t", row.names=F, quote=F)

pdf(paste0(gemma_prefix, ".distribution.pdf"), width = 3, height = 3)
ggplot(pip_df, aes(x=PIP)) + geom_histogram(bins=100) + xlab("PIP") + ylab("SNPs")
dev.off()


lim = round(max(clean_pip_df$PIP)+0.1, 1)
top = round(max(clean_pip_df$PIP)-0.2, 1)
pdf(paste0(gemma_prefix,".manhattan_pip.pdf"), width = 7, height = 3)
manhattan(clean_pip_df, chr="chr", bp="ps", snp="rs", p="PIP", suggestiveline=F, genomewideline=F, logp=F, ylab="PIP", ylim=c(0, lim), annotatePval=top)
dev.off()

#hyp_df = read.csv(paste0(gemma_prefix, ".hyp.txt"), sep="\t")
#high_pve_iter = subset(hyp_df, pve > 0.9)

#snp_lists = list()
#for (i in 1:nrow(high_pve_iter)) {
#   iter = rownames(high_pve_iter)[i]
#   nsnp = high_pve_iter[i, "n_gamma"]
#   snp_lists[[iter]] = as.character(unlist(gamma[iter, 1:nsnp]))
#}
#param[snp_lists[[1]],]
