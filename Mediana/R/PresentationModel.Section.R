######################################################################################################################

# Function: PresentationModel.Section
# Argument: Section object.
# Description: This function is called by default if the class of the argument is a Section object.
#' @export
PresentationModel.Section = function(section, ...) {
  presentationmodel = PresentationModel()
  presentationmodel = presentationmodel + section

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      presentationmodel = presentationmodel + args[[i]]
    }
  }
  return(presentationmodel)
}