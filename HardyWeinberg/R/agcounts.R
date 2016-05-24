agcounts <- function(x,verbose=FALSE){
  lab <- names(x)
  if(!all(lab %in% c("A","AA","AB","B","BB")))
    stop("Unknown genotypes occurred. Supply counts as a named vector c(A,AA,AB,B,BB)")
  n <- sum(x)         
  fAA <- x[lab=="AA"]
  fAB <- x[lab=="AB"]
  fBB <- x[lab=="BB"]
  mA <- x[lab=="A"]
  mB <- x[lab=="B"]
  nAf <- 2*fAA + fAB
  nBf <- 2*fBB + fAB
  nm <- mA + mB
  nf <- n - nm
  nA <- mA + 2*fAA + fAB
  nB <- mB + 2*fBB + fAB
  nt <- nA + nB
  pA <- nA/nt
  pB <- 1 - pA
  if(nA > nB) { #swap counts if A is not the minor allele
      nA <- nB
      nB <- nt - nA
      mA <- mB
      mB <- nm - mA
      fAA <- fBB
      fBB <- nf - (fAA + fAB)
      pA <- pB
      pB <- 1 - pA
  }
  if(verbose) {
     cat("n =",n,"nf = ",nf,"nm =",nm,"\n")
     cat("nt =",nt,"nA =",nA,"nB =",nB,"\n")
     cat("nAf =",nAf,"nBf =",nBf,"\n")
     cat("mA =",mA,"mB = ",mB,"\n")
     cat("fAA =",fAA,"fAB =",fAB,"fBB =",fBB,"\n")
     cat("pA =",round(pA,digits=4),"pB =",round(pB,digits=4),"\n")
  }
  out <- list(nA=nA,nB=nB,nf=nf,nm=nm,n=n,nAf=nAf,nBf=nBf,nt=nt,
       fAA=fAA,fAB=fAB,fBB=fBB,mA=mA,mB=mB,pA=pA,pB=pB)
}
