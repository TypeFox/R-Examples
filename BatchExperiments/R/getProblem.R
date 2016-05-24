#' @title Get problem from registry by id.
#'
#' @description
#' The requested object is loaded from disk.
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param id [\code{character(1)}]\cr
#'   Id of problem.
#' @return [\code{\link{Problem}}].
#' @family get
#' @export
getProblem = function(reg, id) {
  checkExperimentRegistry(reg, strict = TRUE)
  assertString(id)
  pids = dbGetAllProblemIds(reg)
  if (id %nin% pids)
    stop("Unknown problem id, possible candidates are: ", collapse(pids))
  loadProblem(reg, id)
}
