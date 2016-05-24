#' @title ExperimentJob
#'
#' @description
#' You can access job properties using the \code{job} object which is optionally passed
#' to dynamic problem functions and algorithms. The object is a named list with the following
#' elements:
#' \describe{
#'   \item{\code{id} [\code{integer(1)}]:}{Job ID.}
#'   \item{\code{prob.id} [\code{character(1)}]:}{Problem ID.}
#'   \item{\code{prob.pars} [\code{list}]:}{Problem parameters as named list.}
#'   \item{\code{algo.id} [\code{character(1)}]:}{algo.id}{Algorithm ID.}
#'   \item{\code{algo.pars} [\code{list}]:}{Algorithm parameters as named list.}
#'   \item{\code{repl} [\code{integer(1)}]:}{Replication number of this experiment.}
#'   \item{\code{seed} [\code{integer(1)}]:}{Seed set right before algorithm execution.}
#'   \item{\code{prob.seed} [\code{integer(1)}]:}{Seed set right before generation of problem instance.}
#' }
#' @name ExperimentJob
#' @rdname ExperimentJob
NULL

makeExperimentJob = function(id = NA_integer_, prob.id, prob.pars, algo.id, algo.pars, repl, seed, prob.seed) {
  setClasses(list(id = id, prob.id = prob.id, prob.pars = prob.pars, algo.id = algo.id,
                  algo.pars = algo.pars, repl = repl, seed = seed, prob.seed = prob.seed),
             c("ExperimentJob", "Job"))
}

#' @export
print.ExperimentJob = function(x, ...) {
  cat("Experiment:", "\n")
  cat("  Problem:", x$prob.id, "\n")
  cat("  Problem parameters:", convertToShortString(x$prob.pars), "\n")
  cat("  Algorithm:", x$algo.id, "\n")
  cat("  Algorithm parameters:", convertToShortString(x$algo.pars), "\n")
  cat("  Replication:", x$repl, "\n")
  cat("  Seed:", x$seed, "\n")
  cat("  Problem seed:", x$prob.seed, "\n")
}
