######################################################################################################################

# Function: EvaluationModel.Criterion
# Argument: Criterion object.
# Description: This function is called by default if the class of the argument is an Criterion object.
#' @export
EvaluationModel.Criterion = function(criterion, ...) {
  evaluationmodel = EvaluationModel()
  evaluationmodel = evaluationmodel + criterion

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      evaluationmodel = evaluationmodel + args[[i]]
    }
  }
  return(evaluationmodel)
}