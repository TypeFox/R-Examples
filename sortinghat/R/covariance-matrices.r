#' Constructs a p-dimensional intraclass covariance matrix.
#'
#' We define a \eqn{p \times p} intraclass covariance (correlation) matrix to be
#' \eqn{\Sigma = \sigma^2 (1 - \rho) J_p + \rho I_p}, where \eqn{-(p-1)^{-1} <
#' \rho < 1}, \eqn{I_p} is the \eqn{p \times p} identity matrix, and \eqn{J_p}
#' denotes the \eqn{p \times p} matrix of ones.
#'
#' @export
#' @param p the size of the covariance matrix
#' @param rho the intraclass covariance (correlation) constant
#' @param sigma2 the coefficient of the intraclass covariance matrix
#' @return an intraclass covariance matrix of size \eqn{p \times p}
cov_intraclass <- function(p, rho, sigma2 = 1) {
  p <- as.integer(p)
  rho <- as.numeric(rho)
  sigma2 <- as.numeric(sigma2)

  if (rho <= -(p-1)^(-1) || rho >= 1) {
    stop("The value for 'rho' must be between -1 / (p - 1) and 1, exclusively.")
  }
  if (sigma2 <= 0) {
    stop("The value for 'sigma2' must be positive.")
  }
  sigma2 * ((1 - rho) * matrix(1, nrow = p, ncol = p) + rho * diag(p))
}

#' Constructs a p-dimensional covariance matrix with an autocorrelation
#' (autoregressive) structure.
#'
#' This function generates a \eqn{p \times p} autocorrelated covariance matrix
#' with autocorrelation parameter \code{rho}. The variance \code{sigma2} is
#' constant for each feature and defaulted to 1.
#'
#' The autocorrelated covariance matrix is defined as:
#' The \eqn{(i,j)}th entry of the autocorrelated covariance matrix is defined as:
#' \eqn{\rho^{|i - j|}}.
#'
#' The value of \code{rho} must be such that \eqn{|\rho| < 1} to ensure that
#' the covariance matrix is positive definite.
#'
#' @export
#' @param p the size of the covariance matrix
#' @param rho the autocorrelation value
#' @param sigma2 the variance of each feature
#' @return autocorrelated covariance matrix
cov_autocorrelation <- function(p = 100, rho = 0.9, sigma2 = 1) {
  p <- as.integer(p)
  rho <- as.numeric(rho)
  sigma2 <- as.numeric(sigma2)
  
  if (abs(rho) >= 1) {
    stop("The value of 'rho' must be between -1 and 1, exclusively.")
  }
  if (sigma2 <= 0) {
    stop("The value of 'sigma2' must be positive.")
  }
  x <- diag(p)
  sigma2 * rho^abs(row(x) - col(x))
}

#' Generates a p-dimensional block-diagonal covariance matrix with
#' autocorrelated blocks.
#'
#' This function generates a \eqn{p \times p} covariance matrix with
#' autocorrelated blocks. The autocorrelation parameter is \code{rho}. There are
#' \code{num_blocks} blocks each with size, \code{block_size}.  The variance,
#' \code{sigma2}, is constant for each feature and defaulted to 1.
#'
#' The autocorrelated covariance matrix is defined as:
#' \deqn{\Sigma = \Sigma^{(\rho)} \oplus \Sigma^{(-\rho)} \oplus \ldots \oplus
#' \Sigma^{(\rho)},} where \eqn{\oplus} denotes the direct sum and the
#' \eqn{(i,j)}th entry of \eqn{\Sigma^{(\rho)}} is \deqn{\Sigma_{ij}^{(\rho)} =
#' \{ \rho^{|i - j|} \}.}
#'
#' The matrix \eqn{\Sigma^{(\rho)}} is the autocorrelated block discussed above.
#'
#' The value of \code{rho} must be such that \eqn{|\rho| < 1} to ensure that
#' the covariance matrix is positive definite.
#'
#' The size of the resulting matrix is \eqn{p \times p}, where
#' \code{p = num_blocks * block_size}.
#'
#' The block-diagonal covariance matrix with autocorrelated blocks was
#' popularized by Guo et al. (2007) for studying classification of
#' high-dimensional data.
#'
#' @references Guo, Y., Hastie, T., & Tibshirani, R. (2007). "Regularized linear
#' discriminant analysis and its application in microarrays," Biostatistics, 8,
#' 1, 86-100.
#' @export
#' @importFrom bdsmatrix bdsmatrix
#' @param num_blocks the number of blocks in the covariance matrix
#' @param block_size the size of each square block within the covariance matrix
#' @param rho the autocorrelation parameter. Must be less than 1 in absolute
#' value.
#' @param sigma2 the variance of each feature
#' @return autocorrelated covariance matrix
cov_block_autocorrelation <- function(num_blocks, block_size, rho, sigma2 = 1) {
  cov_block <- cov_autocorrelation(p = block_size, rho = rho, sigma2 = sigma2)
  cov_block <- as.vector(cov_block)
  as.matrix(bdsmatrix(blocksize = rep(block_size, num_blocks),
                      blocks = replicate(num_blocks, cov_block)))
}
