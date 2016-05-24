######################################################################################################################

# Function: PresentationModel.
# Argument: ....
# Description: This function is used to call the corresponding function according to the class of the argument.
#' @export
PresentationModel = function(...) {
  UseMethod("PresentationModel")
}