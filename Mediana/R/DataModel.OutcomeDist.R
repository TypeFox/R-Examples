######################################################################################################################

# Function: DataModel.OutcomeDist
# Argument: OutcomeDist object.
# Description: This function is called by default if the class of the argument is an OutcomeDist object.
#' @export
DataModel.OutcomeDist = function(outcome, ...) {
  datamodel = DataModel()
  datamodel = datamodel + outcome

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      datamodel = datamodel + args[[i]]
    }
  }
  return(datamodel)
}