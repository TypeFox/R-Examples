#' Generate delta matrix
#'
#' This function will generate the matrix of deltas, as specified in the paper.
#' See Details section.
#' 
#' @param theta.diff A vector of length p-1, where p is the number of populations
#' of treatments. Coordinate [i] in theta.diff corresponds to \eqn{\theta_i -
#' \theta_{i+1}}.
#' @param sigma.2 The known variance of the error terms.
#' @param n The number of replications in each population.
#'
#' @export
#'
#' @details As specified in the paper, we can assume that the thetas are in a
#' decreasing order, meaning that \eqn{\theta_1 \ge \theta_2, \ldots, \theta_n}.
#' It follows that all the components of the theta.diff vector must be positive.
#' Note that the delta matrix in the paper is a scaled version of the
#' differences between the thetas.
#'
#' @examples
#' del1 <- c(2, 4)
#' genDelMat(del1)
#'
#' @seealso
#' \code{\link{exactCoverageProb}}, \code{\link{integrand}}
#'
#' @return The function returns a matrix with p rows and p columns, that
#' contains the \eqn{delta_{ij}}'s, as described in the paper.

genDelMat <- function(theta.diff, sigma.2=1, n=1) {

  # store the lengths of theta.diff and p
  len.d <- length(theta.diff)
  len.p <- len.d + 1
  
  # set up the index matrix, on which apply function will be used.
  # The index matrix for p=3 would be:
  # 1 1 --> sum theta.diff[1] to theta.diff[1]
  # 1 2 --> sum theta.diff[1] to theta.diff[2]
  # 1 3
  # 2 2 
  # 2 3 --> sum theta.diff[2] to theta.diff[3]
  # 3 3
  # The 'apply' function then sums the corresponding theta.diff coordinates.
  a1 <- matrix(rep(1:len.d, times=len.d), nrow=len.d)
  indexMat <- cbind(rep(1:len.d, times=len.d:1), 
    a1[lower.tri(a1, TRUE)])

  orgMatVals <- apply(indexMat, 1, function(x) 
    sum(theta.diff[x[1]:x[2]]))

  # set up the delta matrix
  orgMat <- matrix(0, len.p, len.p)
  orgMat[lower.tri(orgMat)] <- orgMatVals
  out <- -(orgMat - t(orgMat))
  out*sqrt(n)/sqrt(sigma.2)
}
