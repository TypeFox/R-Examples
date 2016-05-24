binom.blaker.limits <- function(x,n,level=.95,tol=1e-10,...) {
  if (n < 1 || x < 0 || x > n) stop("Parameters n = ",n,", x = ",x, " wrong!")
  if (level <= 0 || level >= 1) stop("Confidence level ",level," out of (0, 1)!")
  if (tol <= 0) stop("Numerical tolerance ",tol," nonpositive!")
  lower <- binom.blaker.lower.limit(x,n,level,tol,...)
  upper <- 1 - binom.blaker.lower.limit(n-x,n,level,tol,...)
  return(c(lower,upper))
}
