######################################################################################################################

# Function: PresentationModel.Subsection
# Argument: Subsection object.
# Description: This function is called by default if the class of the argument is a Subsection object.
#' @export
PresentationModel.Subsection = function(subsection, ...) {
  presentationmodel = PresentationModel()
  presentationmodel = presentationmodel + subsection

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      presentationmodel = presentationmodel + args[[i]]
    }
  }
  return(presentationmodel)
}