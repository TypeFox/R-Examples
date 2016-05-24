##
## Methods for returning the orthogonal matrix "M"
## ===============================================
##
## Method definition for objects of class "Orthom"
##
setMethod(f = "M", signature(object = "Orthom"), function(object) object@M)
