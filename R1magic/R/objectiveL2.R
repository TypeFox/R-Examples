#' 
#' Objective function for Tikhinov L2 penalty
#'
#'@author Mehmet Suzen
#'@param x, unknown vector
#'@param T, transform bases
#'@param phi, measurement matrix
#'@param y, measurement vector
#'@param lambda, penalty term
#'@note Thank you Jason Xu of Washington University for pointing out complex number handling
#'
objectiveL2 <- function(x,T,phi,y,lambda) {
  X  <- T %*% x
  # OF <- norm( phi %*% X - y , c("F") ) ^ 2 + lambda * norm( X, c("F")) ^ 2
  OF <-  sum(Mod(phi %*%  X - y)^2) + lambda * sum(Mod(x))^(0.5)
 return (OF)
}
