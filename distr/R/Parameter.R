################################
##
## Class: Parameter
##
################################


## Access Methoden
setMethod("name", "Parameter", function(object) object@name)
## Replace Methoden
setReplaceMethod("name", "Parameter", 
                 function(object, value){ object@name <- value; object})
