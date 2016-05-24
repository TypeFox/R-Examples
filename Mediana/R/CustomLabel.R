######################################################################################################################

# Function: CustomLabel.
# Argument: by.
# Description: This function is used to create an object of class CustomLabel.
#' @export
CustomLabel = function(param, label) {

  # Error checks
  if (!is.character(param)) stop("CustomLabel: param must be character.")
  if (!(param %in% c("sample.size", "event", "outcome.parameter", "design.parameter", "multiplicity.adjustment"))) stop("CustomLabel: param is invalid.")
  if (!is.character(label)) stop("CustomLabel: label must be character.")

  custom.label = list(param = param, label = label)

  class(custom.label) = "CustomLabel"
  return(custom.label)
  invisible(custom.label)
}