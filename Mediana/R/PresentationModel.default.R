######################################################################################################################

# Function: PresentationModel.default
# Argument: Multiple objects.
# Description: This function is called by default.
#' @export
PresentationModel.default = function(...) {
  args = list(...)
  if (length(args) > 0) {
    stop("Presentation Model doesn't know how to deal with the parameters")
  }
  else {

    presentationmodel = structure(
      list(project = list(username = "[Unknown User]", title = "[Unknown title]", description = "[No description]"),
           section.by = NULL,
           subsection.by = NULL,
           table.by = NULL,
           custom.label = NULL),
      class = "PresentationModel")
  }
  return(presentationmodel)
}