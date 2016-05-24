# This is a hidden function of the l2boost package.

# l2boost internal method to get hybrid solution (step size, critical value)

# @param rho.m vector of gradient corellations to this point (m)
# @param corr.x correlation matrix
# @param lr current step direction

get.hybrid.solution <- function(rho.m, corr.x, lr) {
  
  # get RHS of theorem 2
  nu.r <- nu.limit(rho.m, corr.x[[lr]], lr)
  
  # Find the vector of next critical values
  lr.PlusOne <- which.min.ind(nu.r, lr)
  
  # select the maximal value.
  nu.r.max <- max(nu.r[lr.PlusOne], na.rm = TRUE)
  
  return(list(nu.r = nu.r, lr.PlusOne = lr.PlusOne, nu.r.max = nu.r.max))
}
