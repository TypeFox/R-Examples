# Author: Robert J. Hijmans
# contact: r.hijmans@gmail.com
# Date : Febrary 2010
# Version 0.1
# Licence GPL v3



if (!isGeneric("predict")) {
	setGeneric("predict", function(object, ...)
		standardGeneric("predict"))
}	

setMethod('predict', signature(object='ConvexHull'), 
	function(object, x, ext=NULL, mask=FALSE, filename='',  ...) {
	
		nc <- nrow(object@polygons@data)
		if ( extends(class(x), 'Raster'))  {
			if (! mask) {
				x <- raster(x)
			}
			if (! is.null(ext)) { 
				x <- crop(x, ext) 
			}
			
			xx <- rasterize(object@polygons, raster(x), field=1, fun='sum', background=0, mask=FALSE, update=FALSE, updateValue="NA", getCover=FALSE, silent=TRUE, ...)
			if (mask) {
				xx <- mask(xx, x, ...)
			}
			
			fun <- function(x){x / nc }
			xx <- calc(xx, fun=fun, filename=filename, ...)
			return(xx)
			
		} else {
		
			if (! inherits(x, 'SpatialPoints') )  {
				x <- data.frame(x[,1:2])
				colnames(x) <- c('x', 'y')
				coordinates(x) <- ~ x + y
			}
			v <- .pointsInPolygons(x, object@polygons, sum)
			return(v)
			
		}
	}
)

