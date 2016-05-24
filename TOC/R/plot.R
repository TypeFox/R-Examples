# Author: Ali Santacruz
# Date :  March 2015
# Version 1.0
# Licence GPL v3


if (!isGeneric("plot")) {
	setGeneric("plot", function(x, ...)
		standardGeneric("plot"))
}	

setMethod("plot", signature='Roc', 
          definition = function(x, labelThres=FALSE, modelLeg="Model", digits=3, nticks=5, digitsL=1, posL = NULL, offsetL = 0.5, ...)  {
            .plotROC(x, labelThres=labelThres, modelLeg=modelLeg, digits=digits, nticks=nticks, digitsL=digitsL, posL = posL, offsetL = offsetL, ...)
            return(invisible(NULL))
          }
)


setMethod("plot", signature ='Toc', 
          definition = function(x, labelThres=FALSE, modelLeg="Model", digits=3, nticks=5, digitsL=1, posL = NULL, offsetL = 0.5, ...)  {
            .plotTOC(x, labelThres=labelThres, modelLeg=modelLeg, digits=digits, nticks=nticks, digitsL=digitsL, posL = posL, offsetL = offsetL, ...)
            return(invisible(NULL))
          }
)	