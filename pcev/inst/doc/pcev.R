## ----ex1-----------------------------------------------------------------
library(pcev)
set.seed(12345)

Y <- matrix(rnorm(100*20), nrow=100)
X <- rnorm(100)

pcev_out <- computePCEV(Y, X)
pcev_out 

## ----ex2-----------------------------------------------------------------
pcev_out2 <- computePCEV(Y, X, Wilks=TRUE)
pcev_out2

## ----ex3-----------------------------------------------------------------
pcev_out3 <- computePCEV(Y, X, shrink=TRUE)
pcev_out3

## ----data----------------------------------------------------------------
data(methylation)
data(pheno)
data(position)

## ----manPlot, fig.align='center', fig.cap="Manhattan plot"---------------
# Compute univariate p-values
fit <- lm(methylation ~ pheno)
pval <- vapply(summary(fit), function(sum) {
  pvalue <- sum$coef[2,4]
  return(pvalue)
}, numeric(1))

# Manhattan plot univariate
plot(position$Pos/1e6, -log10(pval), xlab="Position (Mb)",
     ylab="-log10 pvalue", pch=19, cex=0.5)
abline(h=-log10(8.3*10^-6), lty=2)

## ----cluster, eval=FALSE-------------------------------------------------
#  # Break the region into sub-regions of 500 kB
#  cl <- bumphunter::clusterMaker(chr=position$Chr,
#                                 pos=position$Pos,
#                                 assumeSorted=TRUE,
#                                 maxGap = 500)
#  
#  # Some blocks are too big... put limit at 30
#  index <- cl
#  maxInd <- max(index) + 1
#  
#  blockLengths <- table(index)
#  while(sum(blockLengths > 30) > 0) {
#  
#    for (j in unique(index)) {
#      p <- length(index[index == j])
#      if (p > 30) {
#        q <- floor(p/2); r <- p - q
#        index[index == j] <- c(rep_len(maxInd, q),
#                               rep_len(maxInd + 1, r))
#        maxInd <- maxInd + 2
#      }
#    }
#    blockLengths <- table(index)
#  }
#  
#  # Re-index so that we have consecutive indices
#  cl <- index
#  index <- cl
#  counter <- 0
#  for(j in sort(unique(cl))) {
#    counter <- counter + 1
#    index[index == j] <- counter
#  }

## ----index---------------------------------------------------------------
data(index)
table(table(index))

## ----output, cache=TRUE, tidy=TRUE-------------------------------------------------------------------
pcev_out <- computePCEV(methylation, covariate = pheno,
                        estimation = "block", 
                        inference = "permutation",
                        index = index, nperm=10)
pcev_out

## ----manPlotVIP, fig.align='center', fig.cap="Manhattan-VIMP plot"-----------------------------------
# Manhattan plot VIMP
BLK_boundaries <- c(11235000, 11385000)
plot(position$Pos/1e6, pcev_out$VIMP, xlab = "Position (Mb)",
     ylab = "Variable Importance", pch = 19, cex = 0.5, 
     ylim = c(0,1))
lines(x = BLK_boundaries/1e6, y = rep_len(0.9,2),
      lwd = 3, col = 'red')

## ----------------------------------------------------------------------------------------------------
sessionInfo()

