######################################################################################################################

# Function: PresentationModel.Table
# Argument: Table object.
# Description: This function is called by default if the class of the argument is a Table object.
#' @export
PresentationModel.Table = function(table, ...) {
  presentationmodel = PresentationModel()
  presentationmodel = presentationmodel + table

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      presentationmodel = presentationmodel + args[[i]]
    }
  }
  return(presentationmodel)
}