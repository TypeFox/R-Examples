######################################################################################################################

# Function: DataModel.
# Argument: ....
# Description: This function is used to call the corresponding function according to the class of the argument.
#' @export
DataModel = function(...) {
  UseMethod("DataModel")
}