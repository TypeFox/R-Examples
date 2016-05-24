info.binomial.kgroup <- function(p, group.size=1) {
###-----------------------------------------------------------------------
###      Returns the fisher information for the raw binomial model
###
### p: For l groups, a vector of length l. The probability of an event for
###     a single trial.  For k groups a vector of length k.
###
### group.size: The relative number of observations in each group
###
### Returns: The information from a single trial which is spread over
###          the several groups.
###-----------------------------------------------------------------------

  if (any(p < 0.000001 | p > .999999)) stop ("p must between 0 and 1")
  ngroups <- length(p)

  lngrpsz <- length(group.size)
  if (lngrpsz == 1) group.size <- rep(1,ngroups) else 
         if (ngroups != lngrpsz)
              stop ("lengths of p and group.size must match")

  hess <- -1/(p*(1-p))

  info <- (group.size * -hess) / sum(group.size)
  if (ngroups > 1) info <- diag(info)

  return(info)
}
