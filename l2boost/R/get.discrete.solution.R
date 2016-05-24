# This is a hidden function of the l2boost package.

# l2boost internal method to get discrete step sized solution (critical value, critical point)
#
# @param rho.m vector of gradient corellations to this point (m)
# @param corr.x correlation matrix
# @param lr current step direction
# @param nu l1 shrinkage parameter
get.discrete.solution <- function(rho.m, corr.x, lr, nu) {
  # get the step size, (multiple steps of size nu)
  M.step <- mstep.long(rho.m, corr.x[[lr]], lr, nu)
  
  # determine the critical point and critical value
  # break-ties using the gradient-correlation
  
  # Critical direction using minimal step distance.
  lr.PlusOne <- which.min.ind(M.step, lr)
  
  # Critical step length (possibly more than one)
  Lr <- unique(M.step[lr.PlusOne])
  
  # break ties to determine specific direction chosen.
  lr.PlusOne <- break.ties(rho.m, corr.x[[lr]], lr, lr.PlusOne, Lr, nu)
  
  return(list(Lr = Lr, lr.PlusOne = lr.PlusOne, M.step = M.step))
}
