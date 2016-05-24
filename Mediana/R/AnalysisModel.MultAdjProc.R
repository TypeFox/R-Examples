######################################################################################################################

# Function: AnalysisModel.MultAdjProc
# Argument: MultAdjProc object.
# Description: This function is called by default if the class of the argument is a MultAdjProc object.
#' @export
AnalysisModel.MultAdjProc = function(multadjproc, ...) {

  analysismodel = AnalysisModel()
  analysismodel = analysismodel + multadjproc

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      analysismodel = analysismodel + args[[i]]
    }
  }
  return(analysismodel)

}