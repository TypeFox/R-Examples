#' Remove Size Variation
#'
#' Removes first principal component effect in a covariance matrix.
#'
#' @param cov.matrix Covariance matrix
#' @details Function sets the first eigen value to zero.
#' @return Altered covariance matrix with no variation on former first principal component
#' @author Diogo Melo, Guilherme Garcia
#' @export
#' @examples
#' cov.matrix <- RandomMatrix(10, 1, 1, 10)
#' no.size.cov.matrix <- RemoveSize(cov.matrix)
#' eigen(cov.matrix)
#' eigen(no.size.cov.matrix)
#' @keywords size
RemoveSize <- function (cov.matrix) {
  cov.matrix.svd  <-  svd(cov.matrix)
  size.eigen.vector <- cov.matrix.svd$u[, 1]
  size.eigen.value <- cov.matrix.svd$d[1]
  size.factor <- size.eigen.vector * sqrt(size.eigen.value)
  cov.matrix.size.removed <- cov.matrix - size.factor %*% t(size.factor)
  return (cov.matrix.size.removed)
}
