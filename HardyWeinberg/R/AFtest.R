AFtest <- function(x,verbose=TRUE,...) {
  lab <- names(x)
  if (!all(lab %in% c("A", "AA", "AB", "B", "BB"))) 
            stop("Unknown genotypes occurred. Supply counts as a named vector c(A,AA,AB,B,BB)")
  nfAA <- x[lab == "AA"]
  nfAB <- x[lab == "AB"]
  nfBB <- x[lab == "BB"]
  nmA <- x[lab == "A"]
  nmB <- x[lab == "B"]
  x <- c(nmA, nmB, nfAA, nfAB, nfBB)
  names(x) <- c("A", "B", "AA", "AB", "BB")      
  mA <- x[1]
  mB <- x[2]
  fA <- 2*x[3] + x[4]
  fB <- 2*x[5] + x[4]
  n <- sum(x)
  nm <- x[1] + x[2]
  nf <- sum(x[3:5])
  nt <- nm + 2*nf
  AC <- matrix(c(mA,mB,fA,fB),ncol=2,byrow=TRUE)
  rownames(AC) <- c("M","F")
  colnames(AC) <- c("A","B")
  if(verbose) {
    cat("Fisher Exact test for equality of allele frequencies for males and females.\n\n")
    cat("Table of allele counts:\n")
    print(AC)
    cat("\n")
  }
  out <- fisher.test(AC,...)
  n <- sum(x)
  if(verbose) {
    cat("Sample of",n,"indivduals with",nt,"alleles.","p-value =",out$p.value)
    
  }
  na <- sum(AC)
  out <- list(AC=AC,pval=out$p.value)
}
