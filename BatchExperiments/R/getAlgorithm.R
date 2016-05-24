#' @title Get algorithm from registry by id.
#'
#' @description
#' The requested object is loaded from disk.
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param id [\code{character(1)}]\cr
#'   Id of algorithm.
#' @return [\code{\link{Algorithm}}].
#' @family get
#' @export
getAlgorithm = function(reg, id) {
  checkExperimentRegistry(reg, strict = TRUE)
  assertString(id, "character")
  aids = dbGetAllAlgorithmIds(reg)
  if (id %nin% aids)
    stop("Unknown algorithm id, possible candidates are: ", collapse(aids))
  loadAlgorithm(reg, id)
}
