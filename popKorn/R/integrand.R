#' Evaluate integrand alone
#'
#' This function will evaluate the integrand in the expression for the exact
#' coverage probability. 
#' 
#' @param z This is a real number between -Inf and Inf.
#' @param mat1 This is the matrix of the delta values. This matrix can be
#' generated using the \code{\link{genDelMat}} function.
#'
#' @export
#'
#' @seealso
#' \code{\link{exactCoverageProb}}, \code{\link{integrate2}}
#'
#' @return The function returns a scalar value.

integrand <- function(z, mat1) {
  tmp <- sapply(z, function(x) pnorm(x-mat1), simplify=FALSE)
  tmp <- lapply(tmp, function(x) apply(x, 1, prod, na.rm=TRUE))
  unlist(lapply(tmp, function(x) sum(x))) * dnorm(z)
}
