######################################################################################################################

# Function: EvaluationModel.
# Argument: ....
# Description: This function is used to call the corresponding function according to the class of the argument.
#' @export
EvaluationModel = function(...) {
  UseMethod("EvaluationModel")
}