######################################################################################################################

# Function: Subsection.
# Argument: by.
# Description: This function is used to create an object of class SubSection.
#' @export
Subsection = function(by) {

  # Error checks
  if (!is.character(by)) stop("Subsection: by must be character.")
  if (!any(by %in% c("sample.size", "event", "outcome.parameter", "design.parameter", "multiplicity.adjustment"))) stop("Subsection: the variables included in by are invalid.")


  subsection.report = list(by = by)

  class(subsection.report) = "Subsection"
  return(subsection.report)
  invisible(subsection.report)
}