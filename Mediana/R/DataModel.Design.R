######################################################################################################################

# Function: DataModel.Design
# Argument: Design object.
# Description: This function is called by default if the class of the argument is a Design object.
#' @export
DataModel.Design = function(design, ...) {
  datamodel = DataModel()
  datamodel = datamodel + design

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      datamodel = datamodel + args[[i]]
    }
  }
  return(datamodel)

}