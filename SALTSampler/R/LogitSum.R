LogitSum <- function(x) {
  #For x=logit(p), this function returns s = log(sum(p)) where the
  #sum of p is less than or equal to 1. Calculations are designed to 
  #preserve accuracy even for values numerically near 0 or 1.
  
  #Check input
  if (!exists("x")) {
    stop("x is not defined")
  }
  if (!is.vector(x)) {
    stop("x is not a vector")
  }
  if (!is.numeric(x)) {
    stop("x is not numeric")
  }
  
  #Take logp and logq for all values
  x <- sort(x, decreasing=TRUE)
  lpq <- LogPq(x)
  
  #Find sum of logp, using different calculation method
  #if the max(x) < 0 or max(x) >= 0 to preserve accuracy
  lpm1 <- lpq$logp[-1]
  if (x[1] < 0){ # avoid calculating exp(lp1) when it can be numerically zero
    lp1 <- lpq$logp[1]
    out <- lp1 + log1p(sum(exp(lpm1 - lp1)))
  } else { # avoid calculating exp(lp1) when it can be numerically one
    lq1 <- lpq$logq[1]
    out <- log1p(-exp(lq1) + sum(exp(lpm1)))
  }
  
  #Return log(sum(p))
  return(out)
}
