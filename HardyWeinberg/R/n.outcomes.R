n.outcomes <- function(x,verbose=FALSE) { # number of possible outcomes under the mA, fAB distribution
  n <- sum(x)
  nm <- x[1] + x[2]
  nf <- n - nm
  nt <- nm + 2*nf
  nA <- x[1] + 2*x[3] + x[4]
  nB <- nt-nA
  fAA <- x[3]
  fAB <- x[4]
  fBB <- x[5]
  mA <- x[1]
  mB <- x[2]
  if(verbose) {
    cat(n,nm,nf,"\n")
    cat(nt,nA,nB,"\n")
    cat(fAA,fAB,fBB,"\n")
    cat(mA,mB,"\n")
  }
  minor.allele.overall <- min(nA,nB)
  major.allele.overall <- max(nA,nB)
  if(verbose) cat("minor.allele.overall ",minor.allele.overall)
  mA.max <- min(minor.allele.overall,nm)
  mB.max <- min(major.allele.overall,nm)
  mA.min <- nm - mB.max
  mB.min <- nm - mA.max
  mA.seq <- seq(mA.min,mA.max) # Increasing sequence for males for the minor allele
  mB.seq <- seq(mB.max,mB.min,-1) # Decreasing sequence for males the major allele

  male.minor.allele.seq.for.females <- minor.allele.overall - mA.seq
  male.major.allele.seq.for.females <- major.allele.overall - mB.seq
  
  minor.allele.females <- pmin(male.minor.allele.seq.for.females,male.major.allele.seq.for.females)
  number.outcomes.for.females <- floor(minor.allele.females/2)+1
  n.out <- sum(number.outcomes.for.females)
  Tab <- cbind(mA.seq,mB.seq,
               male.minor.allele.seq.for.females,
               male.major.allele.seq.for.females,
               number.outcomes.for.females)
  colnames(Tab) <- c("MaleMin","MaleMaj","Fem.malemin","Fem.malemaj","Fem#out")
  if(verbose) {
   print(Tab)
  }
  #    n.out <- sum(floor(pmin(nf,nA-mAseq)/2)+1)
  #    mcat(n.out)
  return(list(n.out=n.out,Tab=Tab))
}

