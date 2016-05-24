######################################################################################################################

# Function: MultAdj.
# Argument: ....
# Description: This function is used to call the corresponding function according to the class of the argument.
#' @export
MultAdj = function(...) {
  UseMethod("MultAdj")
}