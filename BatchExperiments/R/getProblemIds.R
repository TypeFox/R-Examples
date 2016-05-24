#' @title Get ids of problems in registry.
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param ids [code{integer}]\cr
#'   Job ids to restrict returned problem ids to.
#' @return [\code{character}].
#' @family get
#' @export
getProblemIds = function(reg, ids) {
  checkExperimentRegistry(reg, strict = TRUE)
  if (missing(ids))
    return(dbGetAllProblemIds(reg))
  BatchJobs:::checkIds(reg, ids)
  unique(dbGetProblemIds(reg, ids))
}
