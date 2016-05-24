######################################################################################################################

# Function: is.AnalysisModel.
# Argument: an object.
# Description: Return if the object is of class AnalysisModel

is.AnalysisModel = function(arg){
  return(any(class(arg)=="AnalysisModel"))
}