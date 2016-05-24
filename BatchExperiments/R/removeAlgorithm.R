#' @title Remove algorithm from registry.
#'
#' @description
#' THIS DELETES ALL FILES REGARDING THIS ALGORITHM, INCLUDING ALL JOBS AND RESULTS!
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param id [\code{character(1)}]\cr
#'   Id of algorithm.
#' @param force [\code{logical(1)}]\cr
#'   Also remove jobs which seem to be still running.
#'   Default is \code{FALSE}.
#' @return Nothing.
#' @family remove
#' @export
removeAlgorithm = function(reg, id, force = FALSE) {
  checkExperimentRegistry(reg, strict = TRUE)
  BatchJobs:::syncRegistry(reg)
  assertString(id)

  if (id %nin% dbGetAllAlgorithmIds(reg))
    stop("Algorithm not present in registry: ", id)

  info("Removing Experiments from database")
  ids = dbFindExperiments(reg, algo.pattern = id, like = FALSE)
  removeExperiments(reg, ids = ids, force = force)
  info("Removing Algorithm from database")
  dbRemoveAlgorithm(reg, id)

  fn = getAlgorithmFilePath(reg$file.dir, id)
  info("Deleting algorithm file: %s", fn)
  ok = file.remove(fn)
  if (!ok)
    warningf("Could not remove algorithm file: %s", fn)
  invisible(NULL)
}
