######################################################################################################################

# Function: is.EvaluationModel.
# Argument: an object.
# Description: Return if the object is of class EvaluationModel

is.EvaluationModel= function(arg){
  return(any(class(arg)=="EvaluationModel"))
}
