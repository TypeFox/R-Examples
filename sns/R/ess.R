############################################################################# 
# Computes the effective sample size using the algorithm in Section 2.3 of
# the paper (http://arxiv.org/pdf/1011.0175v1.pdf) by Madeline Thompson.
# The algorithm is taken from earlier work on 'Initial Sequence Estimators'
# by multiple authors. 
# 
# Args: 
#   x - matrix object with each sample (possibly multivariate) as a row
# Returns:
#   effective sample sizes for the time series in each col of 'x'
############################################################################# 
ess <- function(x, method = c("coda", "ise")) {
  method <- match.arg(method) # TODO: added per JSS suggestion; must be tested
  
  # recursive call to each column of x
  if (NCOL(x) > 1) {
    return (apply(x, 2, ess, method))
  }
  
  # require minimum length of 3 for array
  M <- length(x)
  if (M < 3) return (NA)

  # return zero for constant arrays
  if (sd(x) == 0.0) return (0.0)
  
  # coda method
  if (method == "coda") return (effectiveSize(x))
  
  # ise method
  autocf <- acf(x, plot = FALSE, lag.max = length(x))$acf
  sum <- 0
  prev <- 0
  for (s in 1:(M-2)) {
    rho <- autocf[s+1]
    if ((prev + rho) <= 0) {
      break
    } else {
      sum <- sum + rho
      prev <- rho
    }
  }
  return (M/(1+2*sum))
}

