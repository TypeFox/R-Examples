#' Compute the distance correlation
#' @description Computes the distance correlation between two random variables.
#' @param X A numeric vector or matrix.
#' @param y A numeric vector or matrix.
#'
#' @details Distance correlation measures the statistical dependence between two random
#' variables or two random vectors. The important property is that the distance correlation
#' is zero if and only if two random variables are independent. The details of computing
#' the distance correlation can be seen in Szekely, Rizzo and Bakirov (2007) or
#' \url{https://en.wikipedia.org/wiki/Distance_correlation}. This function is the same as
#' \code{dcor} in \code{energy} package.
#'
#' @return A sample distance correlation.
#'
#' @author Debin Qiu, Jeongyoun Ahn
#' @references
#' Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007), Measuring and Testing Dependence
#' by Correlation of Distances, \emph{Annals of Statistics}, Vol. 35 No. 6, pp. 2769-2794.
#' @examples
#' X <- matrix(rnorm(200),ncol = 2)
#' y <- X%*%matrix(c(1.5,4),ncol = 1) + rnorm(100) # univariate y
#' distcor(X,y)
#'
#' X <- iris[1:50, 1:4]
#' Y <- iris[51:100, 1:4]   # multiple response y
#' distcor(X,Y)
#' @importFrom stats dist
#' @export
distcor <- function(X,y) {
  n <- nrow(X)
  distX <- as.matrix(dist(X))
  disty <- as.matrix(dist(y))
  colMX <- matrix(colMeans(distX),n,n)
  colMy <- matrix(colMeans(disty),n,n)
  distX <- distX - colMX - t(colMX) + matrix(mean(colMX),n,n)
  disty <- disty - colMy - t(colMy) + matrix(mean(colMy),n,n)
  dCov <- sum(distX*disty)/n^2
  dVarX <- sum(distX*distX)/n^2
  dVary <- sum(disty*disty)/n^2
  return(sqrt(dCov/sqrt(dVarX*dVary)))
}
