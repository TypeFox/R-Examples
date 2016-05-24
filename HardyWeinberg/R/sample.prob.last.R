sample.prob.last <- function (n, nm, mA, nA, nfAB) 
{
  nf <- n - nm
  nt <- nm + 2 * nf
  nB <- nt - nA
  mB <- nm - mA
  nfAA <- 0.5 * (nA - mA - nfAB)
  nfBB <- 0.5 * (nB - mB - nfAB)
  log2 <- log(2)
  K <- lgamma(nA+1) + lgamma(nB+1) + lgamma(nf+1) + lgamma(nm+1) - lgamma(mA+1) - lgamma(mB+1) - lgamma(nt+1)
  quotient <- nfAB * log2 - lgamma(nfAA+1) - lgamma(nfAB+1) - lgamma(nfBB+1)
  prob <- exp(K + quotient)
  return(prob)
}
