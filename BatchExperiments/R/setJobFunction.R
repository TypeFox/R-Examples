#' @export
setJobFunction.ExperimentRegistry = function(reg, ids, fun, more.args = list(), reset = TRUE, force = FALSE) {
  stop("setJobFunction not available for BatchExperiments. Please use addProblem or addAlgorithm with overwrite = TRUE")
}
