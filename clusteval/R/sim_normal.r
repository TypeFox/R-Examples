#' Generates random variates from multivariate normal populations with intraclass
#' covariance matrices.
#'
#' We generate \eqn{n_m} observations \eqn{(m = 1, \ldots, M)} from each of
#' \eqn{M} multivariate normal distributions such that the Euclidean distance
#' between each of the means and the origin is equal and scaled by
#' \eqn{\Delta \ge 0}.
#'
#' Let \eqn{\Pi_m} denote the \eqn{m}th population with a \eqn{p}-dimensional
#' multivariate normal distribution, \eqn{N_p(\mu_m, \Sigma_m)} with mean
#' vector \eqn{\mu_m} and covariance matrix \eqn{\Sigma_m}. Also, let \eqn{e_m}
#' be the \eqn{m}th standard basis vector (i.e., the \eqn{m}th element is 1 and
#' the remaining values are 0). Then, we define
#' \deqn{\mu_m = \Delta \sum_{j=1}^{p/M} e_{(p/M)(m-1) + j}.} Note that \code{p}
#' must be divisible by \code{M}. By default, the first 10 dimensions of
#' \eqn{\mu_1} are set to \code{delta} with all remaining dimensions set to 0,
#' the second 10 dimensions of \eqn{\mu_2} are set to \code{delta} with all
#' remaining dimensions set to 0, and so on.
#'
#' Also, we consider intraclass covariance (correlation) matrices such that
#' \eqn{\Sigma_m = \sigma^2 (1 - \rho_m) J_p + \rho_m I_p}, where
#' \eqn{-(p-1)^{-1} < \rho_m < 1}, \eqn{I_p} is the \eqn{p \times p} identity
#' matrix, and \eqn{J_p} denotes the \eqn{p \times p} matrix of ones.
#'
#' By default, we let \eqn{M = 5}, \eqn{\Delta = 0}, and \eqn{\sigma^2 = 1}.
#' Furthermore, we generate 25 observations from each population by default.
#'
#' For \eqn{\Delta = 0} and \eqn{\rho_m = \rho}, \eqn{m = 1, \ldots, M}, the
#' \eqn{M} populations are equal.
#'
#' @param n a vector (of length M) of the sample sizes for each population
#' @param p the dimension of the multivariate normal populations
#' @param rho a vector (of length M) of the intraclass constants for each
#' population
#' @param delta the fixed distance between each population and the origin
#' @param sigma2 the coefficient of each intraclass covariance matrix
#' @param seed seed for random number generation (If NULL, does not set seed)
#' @return named list containing:
#' \describe{
#'   \item{x:}{A matrix whose rows are the observations generated and whose
#'   columns are the \code{p} features (variables)}
#'   \item{y:}{A vector denoting the population from which the observation in
#'   each row was generated.}
#' }
#' @export
#' @examples
#' data_generated <- sim_normal(n = 10 * seq_len(5), seed = 42)
#' dim(data_generated$x)
#' table(data_generated$y)
#'
#' data_generated2 <- sim_normal(p = 10, delta = 2, rho = rep(0.1, 5))
#' table(data_generated2$y)
#' sample_means <- with(data_generated2,
#'                      tapply(seq_along(y), y, function(i) {
#'                             colMeans(x[i,])
#'                      }))
#' (sample_means <- do.call(rbind, sample_means))
sim_normal <- function(n = rep(25, 5), p = 50, rho = rep(0.9, 5), delta = 0,
                       sigma2 = 1, seed = NULL) {
  # The number of populations
  M <- length(n)

  if (delta < 0) {
    stop("The value for 'delta' must be a nonnegative constant.")
  }
  if (sigma2 <= 0) {
    stop("The value for 'sigma2' must be positive.")
  }
  if (length(n) != length(rho)) {
    stop("The length of the vectors 'n' and 'rho' must be equal.")
  }
  if (p %% M != 0) {
    stop("We require that 'p' be divisible by 'M'")
  }
  if(!is.null(seed)) {
    set.seed(seed)
  }

  # A matrix whose rows are the population means.
  means <- lapply(seq.int(M), function(m) {
    mu_m <- matrix(0, nrow = M, ncol = p / M)
    mu_m[m, ] <- delta
    mu_m
  })
  means <- do.call(cbind, means)

  # A list containing each of the covariance matrices.
  cov_list <- lapply(rho, intraclass_cov, p = p)

  # Generates the data in a list of length M.
  # Then, we rbind the data together.
  x <- lapply(seq_len(M), function(m) {
    rmvnorm(n[m], means[m, ], cov_list[[m]])
  })
  x <- do.call(rbind, x)
  y <- do.call(c, sapply(seq_len(M), function(m) {
    rep.int(m, n[m])
  }, simplify = FALSE))
  list(x = x, y = y)
}

#' Construct an intraclass covariance matrix.
#'
#' We define a \eqn{p \times p} intraclass covariance (correlation)
#' matrix to be \eqn{\Sigma_m = \sigma^2 (1 - \rho) J_p + \rho I_p},
#' where \eqn{-(p-1)^{-1} < \rho < 1}, \eqn{I_p} is the
#' \eqn{p \times p} identity matrix, and \eqn{J_p} denotes the
#' \eqn{p \times p} matrix of ones.
#'
#' @param p the dimension of the matrix
#' @param rho the intraclass covariance (correlation) constant
#' @param sigma2 the coefficient of the intraclass covariance matrix
#' @return an intraclass covariance matrix matrix of size p
intraclass_cov <- function(p, rho, sigma2 = 1) {
  if (rho <= -(p-1)^(-1) || rho >= 1) {
    stop("The value for 'rho' must be exclusively between -1 / (p - 1) and 1.")
  }
  if (sigma2 <= 0) {
    stop("The value for 'sigma2' must be positive.")
  }
  sigma2 * ((1 - rho) * matrix(1, nrow = p, ncol = p) + rho * diag(p))
}

