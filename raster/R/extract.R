# Author: Robert J. Hijmans
# Date : October 2010
# Version 1.0
# Licence GPL v3


if (!isGeneric("extract")) {
	setGeneric("extract", function(x, y, ...)
		standardGeneric("extract"))
}	



setMethod('extract', signature(x='Raster', y='vector'), 
function(x, y, ...){ 
	y <- round(y)
	return( .cellValues(x, y, ...) )
})
