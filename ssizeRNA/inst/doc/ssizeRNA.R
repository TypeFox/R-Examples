## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("knitr")
opts_chunk$set(
  message=FALSE)

## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
    BiocStyle::latex()

## ----loadpackage, echo=FALSE------------------------------------------------------------
library(Biobase)
library(edgeR)

## ---------------------------------------------------------------------------------------
library(ssizeRNA)

## ----sim1_size, fig.height=4.5, fig.width=6---------------------------------------------
set.seed(2016)
size1 <- ssizeRNA_single(nGenes = 10000, pi0 = 0.8, m = 200, mu = 10, 
                         disp = 0.1, logfc = log(2), fdr = 0.05, 
                         power = 0.8, maxN = 20)
size1$ssize

## ----sim1_power-------------------------------------------------------------------------
check.power(m = 14, mu = 10, disp = 0.1, logfc = log(2), sims = 10)

## ----sim2_para--------------------------------------------------------------------------
data(hammer.eset)
counts <- exprs(hammer.eset)[, phenoData(hammer.eset)$Time == "2 weeks"]
counts <- counts[rowSums(counts) > 0,]  ## filter zero count genes
trt <- hammer.eset$protocol[which(hammer.eset$Time == "2 weeks")] 

## average read count in control group for each gene
mu <- apply(counts[, trt == "control"], 1, mean)

## dispersion for each gene
d <- DGEList(counts)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
disp <- d$tagwise.dispersion

## ----sim2_size, fig.height=4.5, fig.width=6---------------------------------------------
set.seed(2016)
size2 <- ssizeRNA_vary(nGenes = 10000, pi0 = 0.8, m = 200, mu = mu, 
                       disp = disp, logfc = log(2), fdr = 0.05, 
                       power = 0.8, maxN = 15, replace = FALSE)
size2$ssize

## ----sim3_size, fig.height=4.5, fig.width=6---------------------------------------------
set.seed(2016)
logfc <- function(x){rnorm(x, log(2), 0.5*log(2))}
size3 <- ssizeRNA_vary(nGenes = 10000, pi0 = 0.8, m = 200, mu = mu, 
                       disp = disp, logfc = logfc, fdr = 0.05, 
                       power = 0.8, maxN = 20, replace = FALSE)
size3$ssize

## ----sim3_power-------------------------------------------------------------------------
check.power(m = 14, mu = mu, disp = disp, logfc = logfc, sims = 10, 
            replace = FALSE)

## ----echo=TRUE, result=TRUE-------------------------------------------------------------
sessionInfo()

