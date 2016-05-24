#' Simulate Data According to a Demographic Model
#'
#' This function simulates a model created with \code{\link{coal_model}}.
#' The model can be extended with features, \code{\link{parameter}s} and
#' loci. Read the 'coala-introduction' vignette for detailed instructions
#' on creating and simulating such models.
#'
#' @param object The coalescent model to be simulated
#' @param nsim currently unused
#' @param seed A random seed that is set before simulation.
#' @param ... currently unused
#' @param pars Values for parameters specified with \code{\link{par_named}} or
#'   \code{\link{par_range}}. Should be a named numeric vector.
#' @param cores The number of cores that the independent repetitions from
#'   \code{nsim} will be distributed on.
#'   Must be \code{1} on Windows, and should also be \code{1} when using R
#'   in a graphical environment (e.g. Rstudio).
#' @return A list of summary statistics.
#' @export
#' @importFrom stats simulate
#' @importFrom assertthat assert_that
#' @importFrom parallel mclapply
#'
#' @examples
#' model <- coal_model(10, 3) +
#'   feat_mutation(5) +
#'   sumstat_sfs() +
#'   sumstat_tajimas_d()
#' simulate(model, nsim = 2)
#'
#' model <- coal_model(c(5,10), 20) +
#'   feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
#'   feat_mutation(par_range('theta', 1, 10)) +
#'   sumstat_jsfs()
#' simulate(model, pars=c(tau = 1, theta = 5))
simulate.coalmodel <- function(object, nsim = 1, seed, ...,
                               pars = numeric(0), cores = 1) {

  if (!missing(seed)) set.seed(seed)
  simprog <- select_simprog(object)
  if (is.null(simprog)) stop("No simulator found")

  results <- mclapply(seq(len = nsim), function(i) {
    current_pars <- prepare_pars(pars, object)
    simprog$simulate(object, current_pars)
  }, mc.set.seed = TRUE, mc.cores = cores)

  if (nsim == 1) return(results[[1]])
  results
}
