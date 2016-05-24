#' Convert Simulation Results to abc's Parameter Format
#'
#' This function creates an object compatible with the \code{param}
#' argument of the \code{\link[abc]{abc}} function from coala's simulation
#' results.
#'
#' @param sim_results The simulation results as returned from
#'        \code{\link[=simulate.coalmodel]{simulate}}.
#' @param model The model used for the simulations.
#'
#' @return A data.frame that can be used as \code{param}
#'         argument of \code{\link[abc]{abc}}.
#' @export
#' @seealso For generating abc's summary statistics format:
#'          \code{\link{create_abc_sumstat}}
#'
#' @examples
#' model <- coal_model(10, 1) +
#'   feat_mutation(par_prior("theta", rnorm(1, 5, .5))) +
#'   sumstat_sfs()
#' sim_results <- simulate(model, nsim = 2)
#' abc_param <- create_abc_param(sim_results, model)
#' print(abc_param)
create_abc_param <- function(sim_results, model) {
  # First take care for special case when nsim was 1
  if (!is.null(sim_results$pars)) {
    sim_results <- list(sim_results)
  }

  # Now create parmeter data.frame
  do.call(rbind, lapply(sim_results, function(x) {
    data.frame(t(x$pars))
  }))
}


#' Convert Simulation Results to abc's Summary Statistic Format
#'
#' This function creates an object compatible with the \code{sumstat}
#' argument of the \code{\link[abc]{abc}} function from coala's simulation
#' results. It converts all summary statistics that are in the simulation
#' results and expects that each of them is a numeric vector.
#' Use transformation functions to convert none vector-valued  statistics
#' (e.g. \code{\link{sumstat_jsfs}}, \code{\link{sumstat_omega}} or
#' \code{\link{sumstat_trees}}) into a vector.
#'
#' @inheritParams create_abc_param
#'
#' @return A data.frame that can be used as \code{sumstat}
#'         argument of \code{\link[abc]{abc}}.
#' @export
#' @seealso For generating abc's parameter format:
#'          \code{\link{create_abc_param}}
#'
#' @examples
#' # Using the SFS:
#' model <- coal_model(10, 1) +
#'   feat_mutation(par_prior("theta", rnorm(1, 5, .5))) +
#'   sumstat_sfs()
#' sim_results <- simulate(model, nsim = 2)
#' abc_sumstat <- create_abc_sumstat(sim_results, model)
#' print(abc_sumstat)
#'
#' # Using the JSFS and converting it into a vector:
#' model <- coal_model(c(10, 10), 1) +
#'   feat_mutation(par_prior("theta", rnorm(1, 5, .5))) +
#'   feat_migration(par_prior("m", rnorm(1, .5, .1)), symmetri = TRUE) +
#'   sumstat_jsfs(transformation = function(jsfs) {
#'     c(sum(jsfs[1, ]), sum(jsfs[, 1]), sum(jsfs[-1, -1]))
#'   })
#' sim_results <- simulate(model, nsim = 2)
#' abc_sumstat <- create_abc_sumstat(sim_results, model)
#' print(abc_sumstat)
create_abc_sumstat <- function(sim_results, model) {
  # First take care for special case when nsim was 1
  if (!is.null(sim_results$pars)) {
    sim_results <- list(sim_results)
  }

  sumstats <- get_summary_statistics(model)
  do.call(rbind, lapply(sim_results, function(x) {
    stats_combined <- do.call(c, lapply(sumstats, function(stat) {
      value <- x[[stat$get_name()]]
      assert_that(is.numeric(value))
      if (!is.vector(value)) {
        warning("Value of the ", stat$get_name(), " statistic is not a vector. ",
                "Converting it to vector, but this may not be correct. ",
                "Better use transformation functions to convert the ",
                "statistics to vectors.")
      }
      value
    }))
  }))
}
