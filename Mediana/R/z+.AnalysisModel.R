######################################################################################################################

# Function: +.AnalysisModel.
# Argument: Two objects (AnalysisModel and another object).
# Description: This function is used to add objects to the AnalysisModel object
#' @export
"+.AnalysisModel" = function(analysismodel, object) {

  if (is.null(object))
    return(analysismodel)
  else if (class(object) == "Test"){
    ntest = length(analysismodel$tests)
    analysismodel$tests[[ntest+1]] = unclass(object)
  }
  else if (class(object) == "Statistic"){
    nstatistic = length(analysismodel$statistics)
    analysismodel$statistics[[nstatistic+1]] = unclass(object)
  }
  else if (class(object) == "Interim"){
    analysismodel$general$interim$interim.analysis = unclass(object)
  }
  else if (class(object) == "MultAdjProc"){
    nmultadj = length(analysismodel$general$mult.adjust)
    analysismodel$general$mult.adjust[[nmultadj + 1]] = list(unclass(object))
  }
  else if (class(object) == "MultAdjStrategy"){
    nmultadj = length(analysismodel$general$mult.adjust)
    analysismodel$general$mult.adjust[[nmultadj + 1]] = list(unclass(object))
  }
  else if (class(object) == "MultAdj"){
    nmultadj = length(analysismodel$general$mult.adjust)
    if (length(object)>1) analysismodel$general$mult.adjust = c(analysismodel$general$mult.adjust, unclass(object))
    else analysismodel$general$mult.adjust[[nmultadj + 1]] = unclass(object)
  }
  else stop(paste0("Analysis Model: Impossible to add the object of class ",class(object)," to the Analysis Model"))

  return(analysismodel)

}
