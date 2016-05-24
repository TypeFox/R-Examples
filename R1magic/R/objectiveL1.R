#' 
#' Objective function for ridge L1 penalty
#'
#'@author Mehmet Suzen
#'@param x, unknown vector
#'@param T, transform bases
#'@param phi, measurement matrix
#'@param y, measurement vector
#'@param lambda, penalty term
#'@note Thank you Jason Xu of Washington University for pointing out complex number handling
#'
objectiveL1 <- function(x,T,phi,y,lambda) {
  # Part of R1Magic by mehmet.suzen@physics.org
  X  <- T %*% x
  # OF <- norm( phi %*% X - y , c("F") ) ^ 2 + lambda * norm( X, c("1"))
  OF <- sum(Mod(phi %*%  X - y)^2) + lambda * sum(Mod(x))
 return (OF)
}

