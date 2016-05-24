######################################################################################################################

# Function: Table.
# Argument: by.
# Description: This function is used to create an object of class Table.
#' @export
Table = function(by) {

  # Error checks
  if (!is.character(by)) stop("Table: by must be character.")
  if (!any(by %in% c("sample.size", "event", "outcome.parameter", "design.parameter", "multiplicity.adjustment"))) stop("Table: the variables included in by are invalid.")

  table.report = list(by = by)

  class(table.report) = "Table"
  return(table.report)
  invisible(table.report)
}