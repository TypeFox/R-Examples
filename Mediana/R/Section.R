######################################################################################################################

# Function: Section.
# Argument: by.
# Description: This function is used to create an object of class Section.
#' @export
Section = function(by) {

  # Error checks
  if (!is.character(by)) stop("Section: by must be character.")
  if (!any(by %in% c("sample.size", "event", "outcome.parameter", "design.parameter", "multiplicity.adjustment"))) stop("Section: the variables included in by are invalid.")

  section.report = list(by = by)

  class(section.report) = "Section"
  return(section.report)
  invisible(section.report)
}