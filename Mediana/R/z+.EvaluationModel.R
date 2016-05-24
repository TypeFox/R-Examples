######################################################################################################################

# Function: +.EvaluationModel.
# Argument: Two objects (EvaluationModel and another object).
# Description: This function is used to add objects to the EvaluationModel object
#' @export
"+.EvaluationModel" = function(evaluationmodel, object) {

  if (is.null(object))
    return(evaluationmodel)
  else if (class(object) == "Criterion"){
    ncriteria = length(evaluationmodel$criteria)
    evaluationmodel$criteria[[ncriteria+1]] = unclass(object)
  }

  else stop(paste0("Evaluation Model: Impossible to add the object of class ",class(object)," to the Evaluation Model"))

  return(evaluationmodel)

}