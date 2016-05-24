#' Generates random variates from K multivariate contaminated normal
#' populations.
#'
#' We generate \eqn{n_k} observations \eqn{(k = 1, \ldots, K)} from each of
#' \eqn{K} multivariate contaminated normal distributions. Let \eqn{N_p(\mu,
#' \Sigma)} denote the p-dimensional multivariate normal distribution with mean
#' vector \eqn{\mu} and positive-definite covariance matrix \eqn{\Sigma}.  Then,
#' let the \eqn{k}th population have a \eqn{p}-dimensional multivariate
#' contaminated normal distribution:
#'
#' \deqn{(1 - \epsilon_k) N_p(\mu_k, \Sigma_k) + \epsilon_k N_p(\mu_k, \kappa_k \Sigma_k),}
#'
#' where \eqn{\epsilon_k \in [0, 1]} is the probability of sampling from a
#' contaminated population (i.e., outlier) and \eqn{\kappa_k \ge 1} determines
#' the amount of scale contamination. The contaminated normal distribution can
#' be viewed as a mixture of two multivariate normal random distributions, where
#' the second has a scaled covariance matrix, which can introduce extreme
#' outliers for sufficiently large \eqn{\kappa_k}.
#' 
#' The number of populations, \code{K}, is determined from the length of the
#' vector of sample sizes, code{n}. The mean vectors and covariance matrices
#' each can be given in a list of length \code{K}. If one covariance matrix is
#' given (as a matrix or a list having 1 element), then all populations share
#' this common covariance matrix. The same logic applies to population means.
#'
#' The contamination probabilities in \code{epsilon} can be given as a numeric
#' vector or a single value, in which case the degrees of freedom is replicated
#' \code{K} times. The same idea applies to the scale contamination in the
#' \code{kappa} argument.
#'
#' By default, \code{epsilon} is a vector of zeros, and \code{kappa} is a vector
#' of ones. Hence, no contamination is applied by default.
#'
#' @param n a vector (of length K) of the sample sizes for each population
#' @param mean a vector or a list (of length K) of mean vectors
#' @param cov a symmetric matrix or a list (of length K) of symmetric covariance
#' matrices.
#' @param epsilon a vector (of length K) indicating the probability of sampling
#' a contaminated population (i.e., outlier) for each population
#' @param kappa a vector (of length K) that determines the amount of scale
#' contamination for each population
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
#' # Generates 10 observations from each of two multivariate contaminated normal
#' # populations with equal covariance matrices. Each population has a
#' # contamination probability of 0.05 and scale contamination of 10.
#' mean_list <- list(c(1, 0), c(0, 1))
#' cov_identity <- diag(2)
#' data <- simdata_contaminated(n = c(10, 10), mean = mean_list,
#'                              cov = cov_identity, epsilon = 0.05, kappa = 10,
#'                              seed = 42)
#' dim(data$x)
#' table(data$y)
#'
#' # Generates 10 observations from each of three multivariate contaminated
#' # normal populations with unequal covariance matrices. The contamination
#' # probabilities and scales differ for each population as well.
#' set.seed(42)
#' mean_list <- list(c(-3, -3), c(0, 0), c(3, 3))
#' cov_list <- list(cov_identity, 2 * cov_identity, 3 * cov_identity)
#' data2 <- simdata_contaminated(n = c(10, 10, 10), mean = mean_list,
#'                               cov = cov_list, epsilon = c(0.05, 0.1, 0.2),
#'                               kappa = c(2, 5, 10))
#' dim(data2$x)
#' table(data2$y)
#' 
simdata_contaminated <- function(n, mean, cov, epsilon = rep(0, K),
                                 kappa = rep(1, K), seed = NULL) {

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

  # If only one contamination probability or contamination scalar is given,
  # replicate K times
  epsilon <- as.numeric(epsilon)
  if (length(epsilon) == 1) {
    epsilon <- rep.int(epsilon, times = K)
  } else if (length(epsilon) != K) {
    stop("The number of 'epsilon' values must equal K.")
  }

  kappa <- as.numeric(kappa)
  if (length(kappa) == 1) {
    kappa <- rep.int(kappa, times = K)
  } else if (length(kappa) != K) {
    stop("The number of 'kappa' values must equal K.")
  }

  # Ensure that the values in epsilon are in [0, 1]
  if (!all(0 <= epsilon & epsilon <= 1)) {
    stop("The contamination probabilities 'epsilon' must be in [0,1].")
  }

  # Ensure that the values in kappa are nonnegative
  if (any(kappa < 0)) {
    stop("The scale contamination values in 'kappa' must be nonnegative.")
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

  # TODO: Move contamination code to its own 'rcontam' function.
  # The code below could be made more efficient. By being in its own function,
  # it would be easier to profile and to improve.

  # Generates the data in a list of length K. Then, rbinds the data together.
  x <- lapply(seq_len(K), function(k) {
    # Draw n[k] Bernoulli random variates
    contam_draws <- rbinom(n[k], size = 1, prob = epsilon[k])

    # Draw n[k] multivariate normal random variates with and without contamination
    normal_uncontam <- rmvnorm(n[k], mean[[k]], cov[[k]])
    normal_contam <- rmvnorm(n[k], mean[[k]], kappa[k] * cov[[k]])

    # Sum uncontaminated and contaminated draws
    (1 - contam_draws) * normal_uncontam + contam_draws * normal_contam 
  })
  x <- do.call(rbind, x)
  
  # Generates class labels for each observation in 'x'
  y <- factor(rep(seq_along(n), n))
  
  list(x = x, y = y)
}

