# Author: Ali Santacruz
# Date :  March 2015
# Version 1.0
# Licence GPL v3


if (!isGeneric("scaling")) {
  setGeneric("scaling", function(x, ...)
    standardGeneric("scaling"))
}	

setMethod("scaling", signature='Toc', 
          def = function(x, scalingFactor, newUnits)  {
            x@table[, 2:3] <- x@table[, 2:3]/scalingFactor
            x@prevalence <- x@prevalence/scalingFactor
            x@population <- x@population/scalingFactor
            x@units <- newUnits                              
            return(x)
          }
)
