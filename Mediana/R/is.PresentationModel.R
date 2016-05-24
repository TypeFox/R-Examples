######################################################################################################################

# Function: is.PresentationModel.
# Argument: an object.
# Description: Return if the object is of class PresentationModel

is.PresentationModel = function(arg){
  return(any(class(arg)=="PresentationModel"))
}