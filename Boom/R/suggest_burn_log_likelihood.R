SuggestBurnLogLikelihood <- function(log.likelihood, fraction = .25) {
  ## Suggests a burn-in period for an MCMC chain based on the log
  ## likelihood values simulated on the last leg of the chain.
  ## Args:
  ##   log.likelihood: The MCMC sample path of log likelihood for a
  ##     model.
  ##   fraction: The fraction of the chain that should be used to
  ##     determine the log likelihood lower bound.
  ##
  ## Returns:
  ##   An iteration number to be used as a burn-in.  This can be 0 if
  ##   'fraction' is zero, which will be the case if no burn-in is
  ##   desired.
  ##
  ## Details:
  ##   Look at the last 'fraction' of the log.likelihood sequence and
  ##   find the minimum value.  Then return the first iteration where
  ##   log.likelihood exceeds this value.
  if (fraction < 0) {
    return(0)
  }
  stopifnot(fraction <= 1.0)
  stopifnot(is.numeric(log.likelihood), length(log.likelihood) > 0)
  FindFirst <- function (logical.sequence) {
    # Returns the position of the first TRUE in a sequence of
    # logicals.
    if (logical.sequence[1]) return(1)
    y <- rle(logical.sequence)
    return(y$lengths[1] + 1)
  }
  cutpoint <- round(fraction * length(log.likelihood))
  stopifnot(cutpoint >= 1)
  min.log.likelihood <- min(tail(log.likelihood, cutpoint))
  burn <- FindFirst(log.likelihood >= min.log.likelihood) - 1
  if (burn < 0) {
    burn <- 0
  }
  return(burn)
}
