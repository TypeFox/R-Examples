
# Calculates the highest density interval for MCMC output.

# This function replaces the call to emp.hpd in TeachingDemos
# Based on Mike's code from way back.

# Not exported, no sanity checks.

hdi <- function(object, credMass=0.95) {
  x <- sort(object)  # also removes NAs
  if(length(x) == 0)
    return(c(lower = NA, upper = NA))
  n <- length(x)
  # exclude <- ceiling(n * (1 - credMass)) # Not always the same as...
  exclude <- n - floor(n * credMass)       # Number of values to exclude
  low.poss <- x[1:exclude]             # Possible lower limits...
  upp.poss <- x[(n - exclude + 1):n]   # ... and corresponding upper limits
  best <- which.min(upp.poss - low.poss)      # Combination giving the narrowest interval
  c(lower = low.poss[best], upper = upp.poss[best])
}

