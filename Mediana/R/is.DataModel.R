######################################################################################################################

# Function: is.DataModel.
# Argument: an object.
# Description: Return if the object is of class DataModel

is.DataModel = function(arg){
  return(any(class(arg)=="DataModel"))
}