# Author: Ali Santacruz
# Date :  April 2015
# Version 1.0
# Licence GPL v3


if (!isGeneric("ROC")) {
	setGeneric("ROC", function(index, boolean, ...)
		standardGeneric("ROC"))
}	

setMethod("ROC", signature=c('numeric', 'numeric'), 
          definition = function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, progress=FALSE)  {
            tocd <- .ROCnosp(index, boolean, mask=mask, nthres=nthres, thres=thres, NAval=NAval, progress=progress)
            return(tocd)
          }
          )

setMethod("ROC", signature=c('RasterLayer', 'RasterLayer'), 
          definition = function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, progress=FALSE)  {
            tocd <- .ROCsp(index, boolean, mask=mask, nthres=nthres, thres=thres, NAval=NAval, progress=progress)
            return(tocd)
          }
          )
	