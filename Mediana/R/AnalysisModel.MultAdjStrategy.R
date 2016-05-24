######################################################################################################################

# Function: AnalysisModel.MultAdjStrategy
# Argument: MultAdjStrategy object.
# Description: This function is called by default if the class of the argument is a MultAdjStrategy object.
#' @export
AnalysisModel.MultAdjStrategy = function(multadjstrategy, ...) {

  analysismodel = AnalysisModel()
  analysismodel = analysismodel + multadjstrategy

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      analysismodel = analysismodel + args[[i]]
    }
  }
  return(analysismodel)

}