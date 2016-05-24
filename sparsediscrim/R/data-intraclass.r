#' Generates data from \code{K} multivariate normal data populations, where each
#' population (class) has an intraclass covariance matrix.
#'
#' This function generates \code{K} multivariate normal data sets, where each
#' class is generated with a constant mean vector and an intraclass covariance
#' matrix. The data are returned as a single matrix \code{x} along with a vector
#' of class labels \code{y} that indicates class membership.
#' 
#' For simplicity, we assume that a class mean vector is constant for each
#' feature. That is, we assume that the mean vector of the \eqn{k}th class is
#' \eqn{c_k * j_p}, where \eqn{j_p} is a \eqn{p \times 1} vector of ones and
#' \eqn{c_k} is a real scalar.
#' 
#' The intraclass covariance matrix for the \eqn{k}th class is defined as:
#' \deqn{\sigma_k^2 * (\rho_k * J_p + (1 - \rho_k) * I_p),}
#' where \eqn{J_p} is the \eqn{p \times p} matrix of ones and \eqn{I_p} is the
#' \eqn{p \times p} identity matrix.
#'
#' By default, with \eqn{sigma_k^2 = 1}, the diagonal elements of the intraclass
#' covariance matrix are all 1, while the off-diagonal elements of the matrix
#' are all \code{rho}.
#' 
#' The values of \code{rho} must be between \eqn{(1 - p)^(-1)} and 1,
#' exclusively, to ensure that the covariance matrix is positive definite.
#'
#' The number of classes \code{K} is determined with lazy evaluation as the
#' length of \code{n}.
#'
#' @importFrom mvtnorm rmvnorm
#' @export
#' @param n vector of the sample sizes of each class. The length of \code{n}
#' determines the number of classes \code{K}.
#' @param p the number of features (variables) in the data
#' @param rho vector of the values of the off-diagonal elements for each
#' intraclass covariance matrix. Must equal the length of \code{n}.
#' @param mu vector containing the mean for each class. Must equal the length of
#' \code{n} (i.e., equal to \code{K}).
#' @param sigma2 vector of variances for each class. Must equal the length of
#' \code{n}. Default is 1 for each class.
#' @return named list with elements:
#' \itemize{
#'   \item \code{x}: matrix of observations with \code{n} rows and \code{p}
#' columns
#'   \item \code{y}: vector of class labels that indicates class membership for
#' each observation (row) in \code{x}.
#' }
#' @examples
#' # Generates data from K = 3 classes.
#' data <- generate_intraclass(n = 3:5, p = 5, rho = seq(.1, .9, length = 3),
#'                             mu = c(0, 3, -2))
#' data$x
#' data$y
#' 
#' # Generates data from K = 4 classes. Notice that we use specify a variance.
#' data <- generate_intraclass(n = 3:6, p = 4, rho = seq(0, .9, length = 4),
#'                             mu = c(0, 3, -2, 6), sigma2 = 1:4)
#' data$x
#' data$y
generate_intraclass <- function(n, p, rho, mu, sigma2 = rep(1, K)) {
  p <- as.integer(p)
  rho <- as.numeric(rho)
  n <- as.integer(n)
  mu <- as.numeric(mu)

  K <- length(n)

  if (length(rho) != K) {
    stop("The length of 'rho' must equal the length of 'n'.")
  } else if(length(mu) != K) {
    stop("The length of 'mu' must equal the length of 'n'.")
  } else if(length(sigma2) != K) {
    stop("The length of 'sigma2' must equal the length of 'n'.")
  }

  x <- lapply(seq_len(K), function(k) {
    rmvnorm(n = n[k], mean = rep(mu[k], p),
            sigma = cov_intraclass(p = p, rho = rho[k], sigma2 = sigma2[k]))
  })
  x <- do.call(rbind, x)
  y <- factor(rep(seq_along(n), n))
  
  list(x = x, y = y)
}

#' Generates a \eqn{p \times p} intraclass covariance matrix
#'
#' This function generates a \eqn{p \times p} intraclass covariance matrix with
#' correlation \code{rho}. The variance \code{sigma2} is constant for each
#' feature and defaulted to 1.
#'
#' The intraclass covariance matrix is defined as:
#' \deqn{sigma2 * (\rho * J_p + (1 - \rho) * I_p),}
#' where \eqn{J_p} is the \eqn{p \times p} matrix of ones and \eqn{I_p} is the
#' \eqn{p \times p} identity matrix.
#'
#' By default, with \code{sigma2 = 1}, the diagonal elements of the intraclass
#' covariance matrix are all 1, while the off-diagonal elements of the matrix
#' are all \code{rho}.
#' 
#' The value of \code{rho} must be between \eqn{(1 - p)^(-1)} and 1,
#' exclusively, to ensure that the covariance matrix is positive definite.
#'
#' @param p the size of the covariance matrix
#' @param rho the value of the off-diagonal elements
#' @param sigma2 the variance of each feature
#' @return intraclass covariance matrix
cov_intraclass <- function(p, rho, sigma2 = 1) {
  p <- as.integer(p)
  rho <- as.numeric(rho)
  
  if (rho <= (1 - p)^(-1) || rho >= 1) {
    stop("The value of 'rho' must be between (1 - p)^(-1) and 1, exclusively.")
  }
  if (sigma2 <= 0) {
    stop("The value of 'sigma2' must be positive.")
  }
  sigma2 * (rho * matrix(1, nrow = p,  ncol = p) + (1 - rho) * diag(p))
}
