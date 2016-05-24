######################################################################################################################

# Function: AnalysisModel.default
# Argument: arguments.
# Description: This function is called by default if the class of the argument is neither a MultAdjust,
#              nor a Interim object.
#' @export
AnalysisModel.default = function(...) {
  args = list(...)
  if (length(args) > 0) {
    stop("Analysis Model doesn't know how to deal with the parameters")
  }
  else {
    analysismodel = structure(list(general = list(interim.analysis = NULL,
                                                  mult.adjust = NULL),
                                   tests = NULL,
                                   statistics = NULL),
                              class = "AnalysisModel")
  }
  return(analysismodel)
}