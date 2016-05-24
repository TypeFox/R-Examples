# This is a hidden function of the l2boost package.

# l2boost internal method to get l2boost-lars limiting solution. This function is detailed in
# section 4.7 of Ehrlinger 2011
# 
# @param rho.m vector of gradient corellations to this point (m)
# @param corr.x correlation matrix
# @param activeS Acrive set
# @param qr.tol the tolerance for detecting linear dependencies in the qr.solve
# @param eps.tol limiting value for dynamic boosting step size
# 
# @references Ehrlinger (2011). Regularization: Stagewise Regression and Bagging, Ph.D. Dissertation, Case Western Reserve University
# 
# @seealso \code{\link{qr.solve}}

get.lars.solution <- function(rho.m, corr.x, activeS, qr.tol, eps.tol) {
  rhom.sign <- sign(rho.m[activeS])
  activeS.sign <- rhom.sign/rhom.sign[1]
  #define the "first" coordinate: calculations are all relative to this
  first.coord <- activeS[1]
  Q.active <- t(sapply(activeS, function(i) {corr.x[[i]][activeS]}))
  gmma <- qr.solve(Q.active, activeS.sign, tol = qr.tol)/activeS.sign
  Dl <- rho.m/rho.m[first.coord]
  Rl.den <- 1
  Rl.num <- rowSums(sapply(1:length(activeS), function(j) {
    corr.x[[activeS[j]]] * activeS.sign[j] * gmma[j]}), na.rm = TRUE)
  Rl <- Rl.num / Rl.den
  num <- Dl - Rl
  #limiting nu for each coordinate relative to the first coordinate
  nu.limit <- (1 - abs(num)/(1 - Rl * sign(num))) / Rl.den
  #break ties randomly: exclude negative nu values!!!
  lr.PlusOne.candidates <- which.min.ind(nu.limit, c(activeS, which(nu.limit < eps.tol)))
  if (length(lr.PlusOne.candidates) > 0) {
    lr.PlusOne <- resample(lr.PlusOne.candidates)
    nu.limit <- nu.limit[lr.PlusOne]
  }
  else {
    #near the end of the path variable can sometimes yield negative limiting nu
    #not sure why this happens: TBD TBD TBD
    lr.PlusOne <- resample(which.min.ind(nu.limit, activeS))
    nu.limit <- eps.tol
  }
  return(list(nu.limit = nu.limit, lr.PlusOne = lr.PlusOne,
              gmma = gmma, active.set.sign = activeS.sign))
}
