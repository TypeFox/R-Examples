getpvals <- function(X,nsim=100,pvaluetype="selome") {
  simpval <- NULL
  for(j in 1:nrow(X)) {
    y <- as.vector(X[j,])
    names(y) <- c("AA","AB","BB")
    n <- sum(y)
    nA <- 2*y[1]+y[2]
    nB <- 2*n - nA
    MaxHet <- min(nA, nB)
    #  mcat(n,nA,nB)
    if(MaxHet%%2 == 0) 
      nAB <- seq(0, MaxHet, 2)
    if(MaxHet%%2 == 1) 
      nAB<- seq(1, MaxHet, 2)
    probs <- HWCondProbAB(n,nA,nAB)$p
    nAA <- (nA-nAB)/2
    nBB <- (nB-nAB)/2
    Xaux <- cbind(nAA,nAB,nBB)
    colnames(Xaux) <- c("AA","AB","BB")
    pvals <- HWExactMat(Xaux,pvaluetype=pvaluetype)$pvalvec
    # sample from the p-value distribution
    simpval <- rbind(simpval,sample(pvals,nsim,replace=TRUE,prob=probs))
  }
  return(simpval)
}
