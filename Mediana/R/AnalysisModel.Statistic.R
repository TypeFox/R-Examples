######################################################################################################################

# Function: AnalysisModel.Statistic
# Argument: Statistic object.
# Description: This function is called by default if the class of the argument is a Statistic object.
#' @export
AnalysisModel.Statistic = function(statistic, ...) {

  analysismodel = AnalysisModel()
  analysismodel = analysismodel + statistic

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      analysismodel = analysismodel + args[[i]]
    }
  }
  return(analysismodel)

}