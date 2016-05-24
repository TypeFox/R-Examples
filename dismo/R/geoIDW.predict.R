# Author: Robert J. Hijmans
# Date : Febrary 2010
# Version 0.2
# Licence GPL v3



if (!isGeneric("predict")) {
	setGeneric("predict", function(object, ...)
		standardGeneric("predict"))
}	

setMethod('predict', signature(object='InvDistWeightModel'), 
	function(object, x, ext=NULL, filename='', mask=FALSE, ...) {
	
		if ( inherits(x, 'Raster'))  {
			if (! mask) {
				x <- raster(x)
			}
			if (! is.null(ext)) { 
				x <- crop(x, ext) 
			}
			if (mask) {
				xx <- interpolate(x, object@model[[1]], debug.level=0, ...)
				xx <- mask(xx, x, filename=filename, ...)
			} else {
				xx <- interpolate(x, object@model[[1]], filename=filename, ...)
			}				
		} else {
			if (! inherits(x, 'SpatialPoints') )  {
				x <- data.frame(x[,1:2])
				colnames(x) <- c('x', 'y')
				coordinates(x) <- ~ x + y
			}
			xx <- predict(object@model[[1]], x, debug.level=0)
			xx <- xx@data[,1]
		}
		return(xx)
	}
)

