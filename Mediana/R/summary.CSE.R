############################################################################################################################

# Function: summary.CSE
# Argument: x, a CSE object.
#' @export
summary.CSE = function(object,...) {
  print(object$simulation.results)
}