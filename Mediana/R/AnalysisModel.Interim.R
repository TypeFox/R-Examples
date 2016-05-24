######################################################################################################################

# Function: AnalysisModel.Interim
# Argument: Interim object.
# Description: This function is called by default if the class of the argument is a Interim object.

AnalysisModel.Interim = function(interim, ...) {

  analysismodel = AnalysisModel()
  analysismodel = analysismodel + interim

  args = list(...)
  if (length(args)>0) {
    for (i in 1:length(args)){
      analysismodel = analysismodel + args[[i]]
    }
  }
  return(analysismodel)

}