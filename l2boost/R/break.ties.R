# This is a hidden function of the l2boost package.

# break ties of critical values in boosting critical direction determination.
# 
# This internal helper function which determines which direction is more favorable when two directions have the same value of rho.m.
# Used by \code{\link{l2boost}}
# 
# @param rho.m gradient correlation of current direction at step m
# @param corr.x vector of gradient correlations for all candidate directions (including the current direction) 
# @param lr index of the current direction 
# @param lr.PlusOne index of the next candidate
# @param Lr length of current direction descent
# @param nu shrinkage parameter value
#
# @return a single coordinate index indicating which coordinate direction for the next move.
# 
break.ties <- function(rho.m, corr.x, lr, lr.PlusOne, Lr, nu) {
  if (length(lr.PlusOne) > 1) {
    rho.m.PlusOne <- rho.m[lr.PlusOne] -
      (1 - (1 - nu)^Lr) * rho.m[lr] * corr.x[lr.PlusOne]
    lr.PlusOne <- lr.PlusOne[resample(which.max.ind(abs(rho.m.PlusOne)))]
  }
  lr.PlusOne
}
