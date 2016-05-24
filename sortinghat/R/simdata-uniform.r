#' Generates random variates from multivariate uniform populations.
#'
#' We generate \eqn{n} observations from each of \eqn{K_0} multivariate uniform
#' distributions such that the Euclidean distance between each of the
#' populations and the origin is equal and scaled by \eqn{\Delta \ge 0}.
#'
#' To define the populations, let \eqn{x = (X_1, \ldots, X_p)'} be a
#' multivariate uniformly distributed random vector such that \eqn{X_j \sim
#' U(a_j^{(k)}, b_j^{(k)})} is an independently distributed uniform random
#' variable with \eqn{a_j^{(k)} < b_j^{(k)}} for \eqn{j = 1, \ldots, p}.
#'
#' For each population, we set the mean of the distribution along one feature to
#' \eqn{\Delta}, while the remaining features have mean 0. The objective is to
#' have unit hypercubes with \eqn{p = K_0} where the population centroids
#' separate from each other in orthogonal directions as \eqn{\Delta} increases,
#' with no overlap for \eqn{\Delta \ge 1}.
#'
#' Hence, let \eqn{(a_k^{k}, b_k^{(k)}) = c(\Delta - 1/2, \Delta + 1/2)}. For
#' the remaining ordered pairs, let \eqn{(a_j^{(k)}, b_j^{(k)}) = (-1/2,
#' 1/2)}.
#'
#' We generate \eqn{n_k} observations from \eqn{k}th population.
#'
#' For \eqn{\Delta = 0}, the \eqn{K_0 = 5} populations are equal.
#'
#' Notice that the support of each population is a unit hypercube with \eqn{p =
#' K_0} features. Moreover, for \eqn{\Delta \ge 1}, the populations are mutually
#' exclusive and entirely separated.
#'
#' @param n a vector (of length \eqn{K_0}) of the sample sizes for each
#' population
#' @param delta the fixed distance between each population and the origin
#' @param seed seed for random number generation. (If \code{NULL}, does not set
#' seed)
#' @return named list containing:
#' \describe{
#'   \item{x:}{A matrix whose rows are the observations generated and whose
#'   columns are the \code{p} features (variables)}
#'   \item{y:}{A vector denoting the population from which the observation in
#'   each row was generated.}
#' }
#' @export
#' @examples
#' data_generated <- simdata_uniform(seed = 42)
#' dim(data_generated$x)
#' table(data_generated$y)
#'
#' data_generated2 <- simdata_uniform(n = 10 * seq_len(5), delta = 1.5)
#' table(data_generated2$y)
#' sample_means <- with(data_generated2,
#'                      tapply(seq_along(y), y, function(i) {
#'                             colMeans(x[i,])
#'                      }))
#' (sample_means <- do.call(rbind, sample_means))
simdata_uniform <- function(n = rep(25, 5), delta = 0, seed = NULL) {
  if (delta < 0) {
    stop("The value for 'delta' must be a nonnegative constant.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- as.integer(n)
  K_0 <- p <- length(n)

  # Helper function that generates 'n' observations from a multivariate uniform
  # distribution with the configuration given in the description above. The
  # 'unif_params' argument is a matrix having 2 rows, such that the jth column
  # contains a_j and b_j as described above.
  multivariate_uniform <- function(n, unif_params) {
    apply(unif_params, 2, function(params) {
      runif(n = n, min = params[1], max = params[2])
    })
  }
  
  x <- lapply(seq_len(K_0), function(k) {
    bounds <- replicate(n = p, c(-1/2, 1/2))
    bounds[, k] <- bounds[, k] + delta
    multivariate_uniform(n[k], bounds)
  })
  x <- do.call(rbind, x)
  y <- do.call(c, sapply(seq_len(K_0), function(k) {
    rep.int(k, n[k])
  }, simplify = FALSE))
  list(x = x, y = y)
}
