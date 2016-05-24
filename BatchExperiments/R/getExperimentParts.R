#' @title Get all parts required to run a single job.
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param id [\code{integer(1)}]\cr
#'   Id of a job.
#' @return [named list]. Returns the \link{Job}, \link{Problem}, \link{Instance} and \link{Algorithm}.
#' @family get
#' @export
getExperimentParts = function(reg, id) {
  checkExperimentRegistry(reg, strict = TRUE)
  id = BatchJobs:::checkId(reg, id)

  res = namedList(c("job", "prob", "instance", "algo"))
  res$job = dbGetJobs(reg, id)[[1L]]
  res$prob = loadProblem(reg, res$job$prob.id)
  # use insert to keep the slot even if this is NULL
  res = insert(res, list(instance = calcDynamic(reg, res$job, res$prob$static, res$prob$dynamic)))
  res$algo = loadAlgorithm(reg, res$job$algo.id)
  return(res)
}
