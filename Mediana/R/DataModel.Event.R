######################################################################################################################

# Function: DataModel.Event
# Argument: Event object.
# Description: This function is called by default if the class of the argument is an Event object.
#' @export
DataModel.Event = function(event, ...) {
  datamodel = DataModel()
  datamodel = datamodel + event

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      datamodel = datamodel + args[[i]]
    }
  }
  return(datamodel)

}