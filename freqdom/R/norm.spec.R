#' Compute the spectral norm of any matrix P, defined as
#' the square root of the largest eigenvalue of \eqn{P P'}.
#'
#' @title Compute a spectral norm of given matrix P
#' @param P numeric matrix
#' @return numeric norm
#' @examples
#' P = matrix(rnorm(15),3,5)
#' norm.spec(P)
#' @import mvtnorm
#' @importFrom stats rnorm
#' @export 
norm.spec = function(P){
  P = matrix(P)
  if (dim (P) > 2 && !is.numeric(P))
    stop("P must be a numeric matrix or scalar")

  sqrt(eigen(P%*%t(P))$values[1])
}
