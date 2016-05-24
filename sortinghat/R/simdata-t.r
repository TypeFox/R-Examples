#' Generates random variates from K multivariate Student's t populations.
#'
#' We generate \eqn{n_k} observations \eqn{(k = 1, \ldots, K_0)} from each of
#' \eqn{K_0} multivariate Student's t distributions such that the Euclidean
#' distance between each of the means and the origin is equal and scaled by
#' \eqn{\Delta \ge 0}.
#'
#' Let \eqn{\Pi_k} denote the \eqn{k}th population with a \eqn{p}-dimensional
#' multivariate Student's t distribution, \eqn{T_p(\mu_k, \Sigma_k, c_k)}, where
#' \eqn{\mu_k} is the population location vector, \eqn{\Sigma_k} is the
#' positive-definite covariance matrix, and \eqn{c_k} is the degrees of freedom.
#'
#' For small values of \eqn{c_k}, the tails are heavier, and, therefore, the
#' average number of outlying observations is increased.
#'
#' The number of populations, \code{K}, is determined from the length of the
#' vector of sample sizes, code{n}. The centroid vectors and covariance matrices
#' each can be given in a list of length \code{K}. If one covariance matrix is
#' given (as a matrix or a list having 1 element), then all populations share
#' this common covariance matrix. The same logic applies to population
#' centroids. The degrees of freedom can be given as a numeric vector or a
#' single value, in which case the degrees of freedom is replicated \code{K}
#' times.
#'
#' @param n a vector (of length K) of the sample sizes for each population
#' @param centroid a vector or a list (of length K) of centroid vectors
#' @param cov a symmetric matrix or a list (of length K) of symmetric covariance
#' matrices.
#' @param df a vector (of length K) of the degrees of freedom for each
#' population
#' @param seed seed for random number generation (If \code{NULL}, does not set
#' seed)
#' \describe{
#'   \item{x:}{A matrix whose rows are the observations generated and whose
#'   columns are the \code{p} features (variables)}
#'   \item{y:}{A vector denoting the population from which the observation in
#'   each row was generated.}
#' }
#' @importFrom mvtnorm rmvt
#' @export
#' @examples
#' # Generates 10 observations from each of two multivariate t populations
#' # with equal covariance matrices and equal degrees of freedom.
#' centroid_list <- list(c(3, 0), c(0, 3))
#' cov_identity <- diag(2)
#' data_generated <- simdata_t(n = c(10, 10), centroid = centroid_list,
#'                             cov = cov_identity, df = 4, seed = 42)
#' dim(data_generated$x)
#' table(data_generated$y)
#'
#' # Generates 10 observations from each of three multivariate t populations
#' # with unequal covariance matrices and unequal degrees of freedom.
#' set.seed(42)
#' centroid_list <- list(c(-3, -3), c(0, 0), c(3, 3))
#' cov_list <- list(cov_identity, 2 * cov_identity, 3 * cov_identity)
#' data_generated2 <- simdata_t(n = c(10, 10, 10), centroid = centroid_list,
#'                              cov = cov_list, df = c(4, 6, 10))
#' dim(data_generated2$x)
#' table(data_generated2$y)
#' 
simdata_t <- function(n, centroid, cov, df, seed = NULL) {
  if (missing(n) || missing(centroid) || missing(cov) || missing(df)) {
    stop("Each of 'n', 'centroid', 'cov', and 'df' must be provided.")
  }

  if (!is.list(mean)) {
    mean <- list(mean)
  }

  if (!is.list(cov)) {
    cov <- list(cov)
  }

  # The number of populations
  K <- length(n)

  # If only one centroid or covariance matrix are given, replicate K times
  if (length(centroid) == 1) {
    centroid <- replicate(K, centroid)
  } else if (length(centroid) != K) {
    stop("The number of 'centroid' vectors must equal K.")
  }
  if (length(cov) == 1) {
    cov <- replicate(K, cov)
  } else if (length(cov) != K) {
    stop("The number of 'cov' vectors must equal K.")
  }

  # If only one degree of freedom is given, replicate K times
  df <- as.numeric(df)
  if (length(df) == 1) {
    df <- rep.int(df, times = K)
  } else if (length(df) != K) {
    stop("The number of 'df' values must equal K.")
  }

  # Ensure the list of centroids are vectors having equal length
  if (!all(sapply(centroid, is.vector)) || !all_equal(sapply(centroid, length))) {
    stop("The 'centroid' vectors must be of the same length.")
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
    rmvt(n[k], delta = centroid[[k]], sigma = cov[[k]], df = df[k])
  })
  x <- do.call(rbind, x)
  
  # Generates class labels for each observation in 'x'
  y <- factor(rep(seq_along(n), n))
  
  list(x = x, y = y)
}

