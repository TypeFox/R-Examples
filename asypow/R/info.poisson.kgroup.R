info.poisson.kgroup <- function(lambda, group.size=1 ) {
###-----------------------------------------------------------------------
###      Returns the information for the raw poisson model
###
### lambda: For k groups, a vector of length k. The mean of the poisson
###           distribution for each group
###
### group.size: The relative number of observations in each group
###
###
### Returns: The information from a single trial spread over the
###          several groups
###-----------------------------------------------------------------------

  if (any(lambda <= 0)) stop("lambda must be positive")
  ngroups <- length(lambda)

  lngrpsz <- length(group.size)
  if (lngrpsz == 1) group.size <- rep(1,ngroups) else 
    if (ngroups != lngrpsz)
      stop ("lengths of lambda and group.size must match")

  hess <- -1 / lambda
  info <- (group.size * -hess) / sum(group.size)
  if (ngroups > 1) info <- diag(info)

  return(info)
}
