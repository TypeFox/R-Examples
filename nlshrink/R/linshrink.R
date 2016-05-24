#' Linear-shrinkage estimator of population eigenvalues.
#'
#' @param X A data matrix.
#' @param k (Optional) Non-negative integer less than \code{ncol(X)}. If \code{k
#'   == 0} (default), \code{X} is assumed to contain 1 class, which will be
#'   centered. If \code{k >= 1}, \code{X} is assumed to contain \code{k}
#'   classes, each of which has already been centered.
#' @return A numeric vector of length \code{ncol(X)}, containing the population
#'   eigenvalue estimates sorted in ascending order.
#' @description \code{linshrink} estimates the population eigenvalues from the
#'   sample eigenvalues by shrinking each sample eigenvalue towards the global
#'   mean based on a shrinkage factor. Details in referenced publications.
#' @references \itemize{ \item Ledoit, O. and Wolf, M. (2004). A
#'   well-conditioned estimator for large-dimensional covariance matrices.
#'   Journal of Multivariate Analysis, 88(2) \item Ledoit, O. and Wolf, M.
#'   (2016). Numerical Implementation of the QuEST function. arXiv:1601.05870
#'   [stat.CO] }
#' @examples
#' linshrink(X = matrix(rnorm(1e4, mean = 5), nrow = 100, ncol = 100)) # 1 class; will be centered
#' linshrink(X = matrix(rnorm(1e4), nrow = 100, ncol = 100), k = 1) # 1 class; no centering
#' @export
linshrink <- function(X, k = 0)
{
  n <- nrow(X)
  p <- ncol(X)
  if (k == 0) {
    X <- X - tcrossprod(rep(1,n), colMeans(X))
    k = 1
  }

  if (n > k) {
    effn <- n-k
  } else stop("k must be strictly less than nrow(X)")
  S <- crossprod(X)/effn
  lambda <- sort(eigen(S,only.values = TRUE)$values)
  lambda[lambda < 0] <- 0
  if (p > effn) lambda [1:(p-effn)] <- 0
  Z <- X*X
  phi <- sum(crossprod(Z) / effn - 2*crossprod(X)* S / effn + S*S)
  gamma <- sum((S - mean(diag(S))*diag(p))^2) #Frobenius norm^2
  shrinkage <- max(0, min(1,phi/gamma/effn))
  lambda_mean <- mean(lambda)
  return (lambda_mean + sqrt(1 - shrinkage)*(lambda - lambda_mean))
}

#' Linear-shrinkage estimator of population covariance matrix.
#'
#' @param X A data matrix.
#' @param k (Optional) Non-negative integer less than \code{ncol(X)}. If \code{k
#'   == 0} (default), \code{X} is assumed to contain 1 class, which will be
#'   centered. If \code{k >= 1}, \code{X} is assumed to contain \code{k}
#'   classes, each of which has already been centered.
#' @return Population covariance matrix estimate. A square positive
#'   semi-definite matrix of dimension \code{ncol(X)}.
#' @description The linear shrinkage estimator of the population covariance
#'   matrix is computed by shrinking the sample covariance matrix towards the
#'   identity matrix based on a shrinkage factor. Note that the eigenvalues of
#'   the population covariance matrix estimate are not the same as the linear
#'   shrinkage estimates of population eigenvalues. Details in referenced
#'   publication.
#' @references \itemize{ \item Ledoit, O. and Wolf, M. (2004). A
#'   well-conditioned estimator for large-dimensional covariance matrices.
#'   Journal of Multivariate Analysis, 88(2)}
#' @examples
#' linshrink_cov(X = matrix(rnorm(1e4, mean = 5), nrow = 100, ncol = 100)) # 1 class; will be centered
#' linshrink_cov(X = matrix(rnorm(1e4), nrow = 100, ncol = 100), k = 1) # 1 class; no centering
#' @export
linshrink_cov <- function(X, k = 0) {
  n <- nrow(X); p <- ncol(X)
  if (k == 0) {
    X <- X - tcrossprod(rep(1,n), colMeans(X))
    k = 1
  }

  if (n > k) effn <- n-k
  else stop("k must be strictly less than nrow(X)")

  S <- crossprod(X)/effn
  Ip <- diag(p)
  m <- sum(S*Ip) / p
  d2 <- sum( (S - m*Ip)^2 ) / p
  b_bar2 <- 1/(p*effn^2) * sum(apply(X, 1, function(x) sum((tcrossprod(x) - S)^2 )))
  b2 <- min(d2, b_bar2)
  a2 <- d2 - b2
  return (b2/d2 * m * Ip + a2 / d2 * S)
}
