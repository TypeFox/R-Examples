######################################################################################################################

# Function: DataModel.default
# Argument: Multiple character strings.
# Description: This function is called by default if the class of the argument is neither an Outcome,
#              nor a SampleSize object.
#' @export
DataModel.default = function(...) {
  args = list(...)
  if (length(args) > 0) {
    stop("Data Model doesn't know how to deal with the parameters")
  }
  else {
    datamodel = structure(list(general = list(outcome.dist = NULL,
                                              outcome.type = NULL,
                                              sample.size = NULL,
                                              event = NULL,
                                              rando.ratio = NULL,
                                              design = NULL),
                               samples = NULL),
                          class = "DataModel")
  }
  return(datamodel)
}