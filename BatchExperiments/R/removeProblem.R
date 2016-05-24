#' @title Remove problem from registry.
#'
#' @description
#' THIS DELETES ALL FILES REGARDING THIS PROBLEM, INCLUDING ALL JOBS AND RESULTS!
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param id [\code{character(1)}]\cr
#'   Id of problem.
#' @param force [\code{logical(1)}]\cr
#'   Also remove jobs which seem to be still running.
#'   Default is \code{FALSE}.
#' @return Nothing.
#' @family remove
#' @export
removeProblem = function(reg, id, force = FALSE) {
  checkExperimentRegistry(reg, strict = TRUE)
  BatchJobs:::syncRegistry(reg)
  assertString(id)

  if (id %nin% dbGetAllProblemIds(reg))
    stop("Problem not present in registry: ", id)

  info("Removing Experiments from database")
  ids = dbFindExperiments(reg, prob.pattern = id, like = FALSE)
  removeExperiments(reg, ids = ids, force = force)
  info("Removing Problem from database")
  dbRemoveProblem(reg, id)

  fn = getProblemFilePaths(reg$file.dir, id)
  info("Deleting problem files: ", collapse(fn, sep = ", "))
  ok = file.remove(fn)
  if (!all(ok))
    warningf("Could not remove problem files: %s", collapse(fn[!ok], sep = ", "))
  invisible(NULL)
}
