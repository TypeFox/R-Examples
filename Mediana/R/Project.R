######################################################################################################################

# Function: Project.
# Argument: username, title, project.
# Description: This function is used to create an object of class Project.
#' @export
Project = function(username = "[Unknown User]", title = "[Unknown title]", description = "[No description]") {

  # Error checks
  if (!is.character(username)) stop("Project: username must be character.")
  if (!is.character(title)) stop("Project: title must be character.")
  if (!is.character(description)) stop("Project: description must be character.")

  project = list(username = username, title = title, description = description)

  class(project) = "Project"
  return(project)
  invisible(project)
}