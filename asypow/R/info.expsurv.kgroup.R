info.expsurv.kgroup <- function(w, L, group.size=1) {
###-----------------------------------------------------------------------
###      Returns the fisher information for the raw exponential
###      survival model
###
### L : The length of the study. If L is length 1 each group has
###     a study of length L. If L is the same length as w
###     each group has a different length.
###
### w : For l groups a vector of length l. The rate, or inverse
###          of the mean of the exponential distribution.
###
### group.size: The relative number of observations in each group
###
### Returns: The information from a single trial which is spread over
###          the several groups.
###-----------------------------------------------------------------------
  ngroups <- length(w)

  if (length(L) == 1) L <- rep(L,ngroups)
  if (length(L) != ngroups)
     stop("lengths of w and L must match")

  lngrpsz <- length(group.size)
  if (lngrpsz == 1) group.size <- rep(1,ngroups) else 
  if (ngroups != lngrpsz)
      stop ("lengths of w and group.size must match")

  nhess <- (exp(-w*L) - 1)/(L*w^3) + 1/(w^2)
  info <- (group.size * nhess) / sum(group.size)
  if (ngroups > 1) info <- diag(info)

  return(info)
}
