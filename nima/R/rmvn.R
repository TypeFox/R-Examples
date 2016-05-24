#' Simulate Multivariate Normal
#'
#' Simulate random variables from a multivariate normal distribution.
#'
#' @param n Number of simulation replicates.
#' @param mu Mean vector.
#' @param sigma Variance-covariance matrix.
#'
#' @details
#' Use Cholesky decomposition of \code{sigma}, from \code{\link[base]{chol}}.
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @return
#' Matrix of size n by \code{length(mu)}, each row corresponding to a replicate.
#'
#' @examples
#' x <- rmvn(100, c(1,2), matrix(c(1,1,1,4), ncol = 2))
#'
#' @seealso
#' \code{\link[stats]{rnorm}}

rmvn <- function(n, mu = 0, sigma = matrix(1)) {
    p <- length(mu)
    if(any(is.na(match(dim(sigma), p)))) {
      stop("Dimension problem!")
    }
    D <- chol(sigma)
    matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p))
}
