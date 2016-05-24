#Source: boot.Cenarma.R
#parametric bootstrap for Cenarma object
#currently only NA case
#
boot.Cenarma <- function(obj, ...) {
  y <- obj$y
  iy <- obj$iy
  n <- length(y)
  ans <- obj$outarima
  betaH <- coef(ans)
  p <- obj$p
  q <- obj$q
  include.mean <- obj$include.mean
  cL <- cR <- rep(NA, n)
  ar <- ma <- numeric(0)
  if(p!=0) ar <- betaH[1:p]
  if(q!=0) ma <-  betaH[(p+1):(p+q)]
  mu <- ifelse(include.mean, betaH[p+q+1], 0)
  sigmaSqAH <- ans$sigma2
  siga <- sqrt(sigmaSqAH)
  Y <- Z <- mu+siga*as.vector(arima.sim(model=list(ar=ar, ma=ma), n=n))
  Y[is.na(y)] <- NA 
  out <- list(y=Y, iy=iy, censorPts=matrix(c(cL,cR), ncol=2), z=Z)
  class(out) <- "cents"
  out
}
  
  
