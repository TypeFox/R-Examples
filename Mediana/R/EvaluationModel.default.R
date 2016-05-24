######################################################################################################################

# Function: EvaluationModel.default
# Argument: Multiple objects.
# Description: This function is called by default if the class of the argument is not a Criterion object.
#' @export
EvaluationModel.default = function(...) {
  args = list(...)
  if (length(args) > 0) {
    stop("Evaluation Model doesn't know how to deal with the parameters")
  }
  else {
    evaluationmodel = structure(list(general = NULL, criteria = NULL),
                                class = "EvaluationModel")
  }
  return(evaluationmodel)
}