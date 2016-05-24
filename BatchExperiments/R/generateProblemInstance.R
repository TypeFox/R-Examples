#' @title Generate dynamic part of problem.
#'
#' @description
#' Calls the dynamic problem function on the static problem part and
#' thereby creates the problem instance.
#' The seeding mechanism is identical to execution on the slave.
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param id [\code{character(1)}]\cr
#'   Id of job.
#' @return Dynamic part of problem.
#' @aliases Instance
#' @export
generateProblemInstance = function(reg, id) {
  checkExperimentRegistry(reg, strict = TRUE)
  job = getJob(reg, id, check.id = TRUE)
  prob = getProblem(reg, job$prob.id)
  calcDynamic(reg, job, prob$static, prob$dynamic)
}
