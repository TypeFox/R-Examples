######################################################################################################################

# Function: PresentationModel.Project
# Argument: Projet object.
# Description: This function is called by default if the class of the argument is a Project object.
#' @export
PresentationModel.Project = function(project, ...) {
  presentationmodel = PresentationModel()
  presentationmodel = presentationmodel + project

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      presentationmodel = presentationmodel + args[[i]]
    }
  }
  return(presentationmodel)
}