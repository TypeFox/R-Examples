vaf <- function(X,hw=FALSE) {
  if (is.vector(X)) {
      n <- sum(X)
      pA <- af(X)
      pAA <- X[1]/n
  }
  else if (is.matrix(X)) {
      n <- rowSums(X)
      pA <- af(X)
      pAA <- X[,1]/n
  }
  if(!hw) y <- (pA+pAA-2*(pA^2))/(2*n) else y <- pA*(1-pA)/(2*n)
  return(y)
}
