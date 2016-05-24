######################################################################################################################

# Function: AnalysisModel.Test
# Argument: Test object.
# Description: This function is called by default if the class of the argument is a Test object.
#' @export
AnalysisModel.Test = function(test, ...) {

  analysismodel = AnalysisModel()
  analysismodel = analysismodel + test

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      analysismodel = analysismodel + args[[i]]
    }
  }
  return(analysismodel)

}