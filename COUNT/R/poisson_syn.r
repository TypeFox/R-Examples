# Generic Synthetic Rate Parameterized Poisson
# Table 6.20, Hilbe, J.M. Negative Binomial Regression. 2nd ed, Cambridge Univ Press
require(MASS)
poisson_syn  <- function(nobs = 50000, off = 0, xv = c(1, -.5,  1))  {
  p <- length(xv) - 1
  X <- cbind(1, matrix(rnorm(nobs * p), ncol = p))
  xb <- X %*% xv
  exb <- exp(xb + off)              
  py <- rpois(nobs, exb)     
  out <- data.frame(cbind(py, X[,-1]))
  names(out) <- c("py", paste("x", 1:p, sep=""))
  return(out)
}


