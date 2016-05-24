######################################################################################################################

# Function: AnalysisModel.MultAdj
# Argument: MultAdj object.
# Description: This function is called by default if the class of the argument is a MultAdj object.
#' @export
AnalysisModel.MultAdj = function(multadj, ...) {

  analysismodel = AnalysisModel()
  analysismodel = analysismodel + multadj

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      analysismodel = analysismodel + args[[i]]
    }
  }
  return(analysismodel)

}