######################################################################################################################

# Function: PresentationModel.CustomLabel
# Argument: CustomLabel object.
# Description: This function is called by default if the class of the argument is a CustomLabel object.
#' @export
PresentationModel.CustomLabel = function(customlabel, ...) {
  presentationmodel = PresentationModel()
  presentationmodel = presentationmodel + customlabel

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      presentationmodel = presentationmodel + args[[i]]
    }
  }
  return(presentationmodel)
}