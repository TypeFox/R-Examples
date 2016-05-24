######################################################################
# Calculate sample size for a binomial parameter in a pooled setting
# based on a confidence interval width specification.
#
# The formula of Worlund & Taylor (1983) is used.
#
# Parameters:
#   p0    - hypothesized upper bound 
#   alpha - an (1-alpha/2)*100% confidence interval is computed
#   d     - half width of the confidence interval
#   k     - number of samples in each pool
# Details:
######################################################################

ciss.pool.wald <- function(pi0, alpha, d, k) {
  n <- (qnorm(1-alpha/2) * (1-pi0) / (d*k))^2 * ((1-pi0)^(-k) - 1)
  return(ceiling(n))
}

