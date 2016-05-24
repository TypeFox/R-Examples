f.QQconf <- function(nGenes=18, nSim=10000, quantiles=c(0.025, 0.975)) {
  nSim <- nSim-1
##  pmat <- matrix(NA, nrow=nSim, ncol=nGenes)
##  for (i in 1:nSim) {
##    pmat[i,] <- sort(-log(runif(nGenes), base=10))
##  }
  pmat <- matrix(-log(runif(nGenes*nSim), base=10), ncol=nGenes)
  pmat <- apply(pmat, 1, sort)
  nSim.first <- 10
  pmat.first <- pmat[,1:nSim.first]
  pmat.sort <- apply(pmat, 1, sort)
  nq <- length(quantiles)
  p.quant <- matrix(NA, nrow=nq, ncol=nGenes)
  for (i in 1:nq) {
    ind <- quantiles[i] * (nSim+1)
    p.quant[i,] <- pmat.sort[ind,]
  }
  rownames(p.quant) <- quantiles
  p.mean <- apply(pmat.sort, 2, mean)
##  cat("p.mean=", p.mean, "\np.mean.true=", p.mean.true, "\n")
##  p.lower <- pmat[]
  list(first=pmat.first, mean=p.mean, quant=p.quant)
}
