auxiliartable <-
function (x, verbose = FALSE) 
{
  n <- sum(x)
  nm <- x[1] + x[2]
  nf <- n - nm
  nt <- nm + 2 * nf
  nA <- x[1] + 2 * x[3] + x[4]
  nB <- nt - nA
  fAA <- x[3]
  fAB <- x[4]
  fBB <- x[5]
  mA <- x[1]
  mB <- x[2]
  if (verbose) {
    cat(n, nm, nf, "\n")
    cat(nt, nA, nB, "\n")
    cat(fAA, fAB, fBB, "\n")
    cat(mA, mB, "\n")
  }
  minor.allele.overall <- min(nA, nB)
  major.allele.overall <- max(nA, nB)
  if (verbose) 
    cat("minor.allele.overall ", minor.allele.overall)
  mA.max <- min(minor.allele.overall, nm)
  mB.max <- min(major.allele.overall, nm)
  mA.min <- nm - mB.max
  mB.min <- nm - mA.max
  mA.seq <- seq(mA.min, mA.max)
  mB.seq <- seq(mB.max, mB.min, -1)
  male.minor.allele.seq.for.females <- minor.allele.overall - 
    mA.seq
  male.major.allele.seq.for.females <- major.allele.overall - 
    mB.seq
  minor.allele.females <- pmin(male.minor.allele.seq.for.females, 
                               male.major.allele.seq.for.females)
  number.outcomes.for.females <- floor(minor.allele.females/2) + 1
  Tab <- cbind(mA.seq, mB.seq, minor.allele.females)
  colnames(Tab) <- c("MaleMin", "MaleMaj","FemaleMin")
  if (verbose) {
    print(Tab)
  }
  return(Tab = Tab)
}
