#' @title Get ids of algorithms in registry.
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param ids [code{integer}]\cr
#'   Job ids to restrict returned algorithm ids to.
#' @return [\code{character}].
#' @family get
#' @export
getAlgorithmIds = function(reg, ids) {
  checkExperimentRegistry(reg, strict = TRUE)
  if (missing(ids))
    return(dbGetAllAlgorithmIds(reg))
  BatchJobs:::checkIds(reg, ids)
  unique(dbGetAlgorithmIds(reg, ids))
}
