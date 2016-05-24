## ----knitr_options, include=FALSE----------------------------------------
library(knitr)
opts_chunk$set(fig.width=7, fig.height=4.5)
options(digits=4, scipen=5)

## ----load_data_hidden, include=FALSE-------------------------------------
library(qtl)
library(lineup)
data(f2cross)
data(expr1)
data(expr2)
data(pmap)
data(genepos)

## ----load_libraries, eval=FALSE------------------------------------------
#  library(qtl)
#  library(lineup)

## ----load_data_shown, eval=FALSE-----------------------------------------
#  data(f2cross)
#  data(expr1)
#  data(expr2)
#  data(pmap)
#  data(genepos)

## ----scale_expr----------------------------------------------------------
expr1 <- expr1/1000
expr2 <- expr2/1000

## ----summary_expr--------------------------------------------------------
nrow(expr1)
nrow(expr2)

## ----find_commond_ind_expr-----------------------------------------------
eid <- findCommonID(expr1, expr2)
length(eid$first)

## ----find_correlated_genes-----------------------------------------------
cor_ee <- corbetw2mat(expr1[eid$first,], expr2[eid$second,], what="paired")

## ----hist_corr_betw_tissues----------------------------------------------
par(mar=c(5,4,1,1))
hist(cor_ee, breaks=seq(-1, 1, len=101), main="", las=1,
     xlab="Correlation in gene expression between tissues")

## ----distee--------------------------------------------------------------
d_ee <- distee(expr1[,abs(cor_ee)>0.9], expr2[,abs(cor_ee)>0.9], d.method="cor")

## ----plot_distee, fig.height=9-------------------------------------------
par(mar=c(5,4,2,1))
plot(d_ee)

## ----count_small_selfself------------------------------------------------
sum(pulldiag(d_ee) < 0.5)

## ----count_large_selfnonself---------------------------------------------
d_ee_nodiag <- omitdiag(d_ee)
sum( !is.na(d_ee_nodiag) & d_ee_nodiag > 0.5)

## ----summary_distee------------------------------------------------------
summary(d_ee)

## ----plot_expr_dup, fig.width=7------------------------------------------
par(mar=c(5,4,1,1))
plot(expr1["48",], expr1["76",],
     xlab="Sample 48, expr1", ylab="Sample 76, expr1",
     las=1, pch=21, bg="slateblue", cex=0.7)

## ----drop_expr1_sample48-------------------------------------------------
expr1["76",] <- colMeans(expr1[c("48", "76"),])
expr1 <- expr1[rownames(expr1)!="48",]

## ----calc_genoprob-------------------------------------------------------
f2cross <- calc.genoprob(f2cross, step=1, error.prob=0.002)

## ----find_pseudomarkers, warning=FALSE-----------------------------------
pmar <- find.gene.pseudomarker(f2cross, pmap, genepos)

## ----calc_locallod, message=FALSE----------------------------------------
id1 <- findCommonID(f2cross, expr1)
lod1 <- calc.locallod(f2cross[,id1$first], expr1[id1$second,], pmar, verbose=FALSE)
id2 <- findCommonID(f2cross, expr2)
lod2 <- calc.locallod(f2cross[,id2$first], expr2[id2$second,], pmar, verbose=FALSE)

## ----disteg--------------------------------------------------------------
d_eg_1 <- disteg(f2cross, expr1[, lod1>25], pmar, verbose=FALSE)
d_eg_2 <- disteg(f2cross, expr2[, lod2>25], pmar, verbose=FALSE)

## ----summary_g_vs_expr1--------------------------------------------------
summary(d_eg_1)

## ----summary_g_vs_expr2--------------------------------------------------
summary(d_eg_2)

## ----combinedist---------------------------------------------------------
d_eg <- combinedist(d_eg_1, d_eg_2)
summary(d_eg)

