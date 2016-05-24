#' Generates random variates from K multivariate normal populations.
#'
#' We generate \eqn{n_k} observations \eqn{(k = 1, \ldots, K)} from each of
#' \eqn{K} multivariate normal distributions. Let the \eqn{k}th population have
#' a \eqn{p}-dimensional multivariate normal distribution, \eqn{N_p(\mu_k,
#' \Sigma_k)} with mean vector \eqn{\mu_k} and positive-definite covariance
#' matrix \eqn{\Sigma_k}.
#'
#' The number of populations, \code{K}, is determined from the length of the
#' vector of sample sizes, code{n}. The mean vectors and covariance matrices
#' each can be given in a list of length \code{K}. If one covariance matrix is
#' given (as a matrix or a list having 1 element), then all populations share
#' this common covariance matrix. The same logic applies to population means.
#'
#' @param n a vector (of length K) of the sample sizes for each population
#' @param mean a vector or a list (of length K) of mean vectors
#' @param cov a symmetric matrix or a list (of length K) of symmetric covariance
#' matrices.
#' @param seed seed for random number generation (If \code{NULL}, does not set
#' seed)
#' @return named list containing:
#' \describe{
#'   \item{x:}{A matrix whose rows are the observations generated and whose
#'   columns are the \code{p} features (variables)}
#'   \item{y:}{A vector denoting the population from which the observation in
#'   each row was generated.}
#' }
#' @importFrom mvtnorm rmvnorm
#' @export
#' @examples
#' 
#' # Generates 10 observations from each of two multivariate normal populations
#' # with equal covariance matrices.
#' mean_list <- list(c(1, 0), c(0, 1))
#' cov_identity <- diag(2)
#' data_generated <- simdata_normal(n = c(10, 10), mean = mean_list,
#'                                  cov = cov_identity, seed = 42)
#' dim(data_generated$x)
#' table(data_generated$y)
#'
#' # Generates 10 observations from each of three multivariate normal
#' # populations with unequal covariance matrices.
#' set.seed(42)
#' mean_list <- list(c(-3, -3), c(0, 0), c(3, 3))
#' cov_list <- list(cov_identity, 2 * cov_identity, 3 * cov_identity)
#' data_generated2 <- simdata_normal(n = c(10, 10, 10), mean = mean_list,
#'                                   cov = cov_list)
#' dim(data_generated2$x)
#' table(data_generated2$y)
#' 
simdata_normal <- function(n, mean, cov, seed = NULL) {
  if (missing(n) || missing(mean) || missing(cov)) {
    stop("Each of 'n', 'mean', and 'cov' must be provided.")
  }

  if (!is.list(mean)) {
    mean <- list(mean)
  }

  if (!is.list(cov)) {
    cov <- list(cov)
  }

  # The number of populations
  K <- length(n)

  # If only one mean or covariance matrix are given, replicate K times
  if (length(mean) == 1) {
    mean <- replicate(K, mean)
  } else if (length(mean) != K) {
    stop("The number of 'mean' vectors must equal K.")
  }
  if (length(cov) == 1) {
    cov <- replicate(K, cov)
  } else if (length(cov) != K) {
    stop("The number of 'cov' vectors must equal K.")
  }

  # Ensure the list of means are vectors having equal length
  if (!all(sapply(mean, is.vector)) || !all_equal(sapply(mean, length))) {
    stop("The 'mean' vectors must be of the same length.")
  }

  # Ensure the list of covariance matrices are square matrices having equal
  # dimensions
  if (!all(sapply(cov, is.matrix)) || !sapply(cov, isSymmetric) ||
      !all_equal(sapply(cov, nrow))) {
    stop("The 'cov' matrices must be symmetric and have the same dimensions.")
  }

  # Sets the RNG seed if provided.
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generates the data in a list of length K. Then, rbinds the data together.
  x <- lapply(seq_len(K), function(k) {
    rmvnorm(n[k], mean[[k]], cov[[k]])
  })
  x <- do.call(rbind, x)
  
  # Generates class labels for each observation in 'x'
  y <- factor(rep(seq_along(n), n))
  
  list(x = x, y = y)
}

