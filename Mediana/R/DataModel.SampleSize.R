######################################################################################################################

# Function: DataModel.SampleSize
# Argument: SampleSize object.
# Description: This function is called by default if the class of the argument is an SampleSize object.
#' @export
DataModel.SampleSize = function(sample.size, ...) {
  datamodel = DataModel()
  datamodel = datamodel + sample.size

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      datamodel = datamodel + args[[i]]
    }
  }
  return(datamodel)
}