# Author: Ali Santacruz
# Date :  April 2015
# Version 1.0
# Licence GPL v3


if (!isGeneric("TOC")) {
	setGeneric("TOC", function(index, boolean, ...)
		standardGeneric("TOC"))
}	

setMethod("TOC", signature=c('numeric', 'numeric'), 
          definition = function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, progress=FALSE, units=character(0)){
            tocd <- .TOCnosp(index, boolean, mask=mask, nthres=nthres, thres=thres, NAval=NAval, P=P, Q=Q, progress=progress, units=units)
            return(tocd)
          }
          )

setMethod("TOC", signature=c('RasterLayer', 'RasterLayer'), 
          definition = function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, progress=FALSE)  {
            tocd <- .TOCsp(index, boolean, mask=mask, nthres=nthres, thres=thres, NAval=NAval, P=P, Q=Q, progress=progress)
            return(tocd)
          }
          )
	