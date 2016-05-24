######################################################################
# Sample size based on mid-p confidence intervals
#
# Params:
#   p0   - hypothesized upper bound (if below 0.5, if above 0.5 then
#          lower bound) on the parameter p.
#   d    - half width of the confidence interval
#  alpha - level of the confidence interval
######################################################################

ciss.midp <- function(p0, d, alpha,nMax=1e6) {
  #Determined by p0 and d
  pi.L <- p0 - d
  pi.U <- p0 + d
  if (pi.L < 0) stop("p0 - d is below zero!")
  if (pi.U > 1) stop("p0 + d is above one!")
  
  #Start from n=1
  n <- floor( max(1/p0, 1/(1-p0))) #0
  done <- FALSE
  #Increase n until the lhs is below alpha
  while (!done & (n<nMax)) {
    n <- n+1
    #Nearest x so x/n \approx p0
    x <- round(p0 * n)
    #LHS
    lhs2 <-  1/2*dbinom(x, size=n, prob=pi.L) +
      1/2*dbinom(x, size=n, prob=pi.U) +
        pbinom(x, size=n, prob=pi.L, lower.tail=FALSE) +
          pbinom(x-1, size=n, prob=pi.U)
    #Are we below alpha already?
    if (!is.na(lhs2)) { done <- (lhs2 < alpha) }
  }
 
  return(n)
}
