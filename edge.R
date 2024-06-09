library(BiocManager)
library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(tidyverse)
setwd("/home/sarath/R/rwd/edge")
count_matrix <- read.csv("matrixcounts.csv", header = T)
head(count_matrix)
View(count_matrix)

sample <- read.csv("sample.csv", header = TRUE)
head(sample)

count_matrix <- count_matrix %>% mutate(across(where(is.numeric), ~replace_na(., median(., na.rm=TRUE))))

dge_list <- DGEList(counts = count_matrix, genes=row.names(count_matrix), group = factor(c("control","treated","control","treated","control","treated","control","treated","control","treated","control","treated")))  
dge_list
keep <- filterByExpr(dge_list)
table(keep)
dge_list <- dge_list[keep, , keep.lib.sizes = F]
dge_list
write.csv(dge_list,"dgelistnew.csv")
dge_list <- calcNormFactors(dge_list)
dge_list$samples
plotMDS(dge_list)

dge_list <- estimateDisp(dge_list)
dge_list$common.dispersion
dge_list$trended.dispersion
dge_list$tagwise.dispersion
dge_list$AveLogCPM
dge_list$trend.method
dge_list$prior.df
dge_list$prior.n
dge_list$span
plotBCV(dge_list)
qlfit <- glmQLFit(dge_list)
qlftest <- glmQLFTest(qlfit)
topTags(qlftest)
write.csv(qlftest,"qlftestnew.csv")
pval <- qlftest[["table"]]
pval <- filter(pval, PValue<0.05)
View(pval)
write.csv(pval,"pvalnew.csv")

fit <-glmFit(dge_list)
lrt <-glmLRT(fit)
topTags(lrt)
view(lrt)
write.csv(lrt,"lrtnew.csv")
pval2 <- lrt[["table"]]
pval2 <- filter(pval2, PValue<0.05)
View(pval2)
write.csv(pval2,"pval_2new.csv")
et <- exactTest(dge_list)
top_degs = topTags(et)
top_degs
write.csv(top_degs,"topdegsct2new.csv")
cpm(dge_list)[rownames(top_degs),]
summary(de <- decideTests(et, lfc = 1))
detags <- rownames(dge_list)[as.logical(de)]
plotSmear(et, de.tags = detags)
abline(h = c(-1,1), col = "blue")

go <- goana(qlftest)
topGO(go,sort = "up")
keg <- kegga(qlftest)
topKEGG(keg, sort = "up")

hist(et$table$PValue, col = "orange")
head(et$table)

### Volcano Plot
plot(et$table$logFC, -10*log10(et$table$PValue), main="Volcano plot", xlab="logFC", ylab="-10*log(P-val)")
Upregulated <- et$table$logFC > 0 & et$table$PValue < 0.01
# highlight our DE genes
points(et$table$logFC[Upregulated], -10*log10(et$table$PValue[Upregulated]), col="red")
# identify genes enriched in the control
downregulated <- et$table$logFC < 0 & et$table$PValue < 0.01
points(et$table$logFC[downregulated], -10*log10(et$table$PValue[downregulated]), col="green")
legend("topright", 
       legend=c("UP", "Down", "Not significant"), 
       col=c("red", "green", "black"), 
       pch=20)

