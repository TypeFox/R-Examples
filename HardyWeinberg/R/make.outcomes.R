make.outcomes <- function(x,verbose=FALSE) {
  #
  # Generates all possible outcomes of the joint (mA,fAB) distribution for the given genotype vector.
  #
  #    x <- c(3,  2,  0,  4,  1) 
  n <- sum(x) # total sample size
  nm <- x[1] + x[2] # number of males
  nf <- n - nm # number of females
  nt <- nm + 2*nf # total numboer of alleles
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
  MA <- rep(mA.seq,Tab[,5])
  NAB <- rep(Tab[,3],Tab[,5]) # repeat female minor.male allele 
  RES <- cbind(MA,nm-MA)
  loose <- pmin(Tab[,3],Tab[,4]) # the minor allele of the females.
#  print(max(loose))
#  print(loose)
  FAB <- NULL
#  fab <- numeric(length(MA))
  for(i in 1:length(loose)) {
#    cat(loose[i], " ")
    if((loose[i]%%2)==0) ss <- seq(0,loose[i],2) # if minor allele even then 0,2
    if((loose[i]%%2)==1) ss <- seq(1,loose[i],2) # if minor allele odd then 1,3,..
#    print(ss)
    FAB <- c(FAB,ss)
  }
#  FAB <- sequ
  FAA <- (NAB-FAB)/2
  FBB <- nf - (FAA+FAB)
  RES <- cbind(RES,FAA,FAB,FBB)
  colnames(RES) <- c("A","B","AA","AB","BB")
  return(RES)
}

