# Author: Robert J. Hijmans
# Date : February 2010
# Version 0.2
# Licence GPL v3


setMethod('predict', signature(object='CirclesRange'), 
	function(object, x, ext=NULL, mask=FALSE, filename='', ...) {
	
		if ( extends(class(x), 'Raster'))  {
			if (! mask) {
				x <- raster(x)
			}
			if (! is.null(ext)) { 
				x <- crop(x, ext) 
			}
			if (mask) {
				xx <- rasterize(object@polygons, raster(x), field=1, fun='sum', mask=FALSE, update=FALSE, updateValue="NA", getCover=FALSE, silent=TRUE, ...)
				xx <- mask(xx, x, filename=filename, ...)
			} else {
				xx <- rasterize(object@polygons, raster(x), field=1, fun='sum', mask=FALSE, update=FALSE, updateValue="NA", getCover=FALSE, silent=TRUE, filename=filename, ...)
			}

			#nc <- length(object@polygons@polygons) 
			#fun <- function(x){x / nc }
			#xx <- calc(xx, fun=fun, filename=filename, progress=progress, ...)
			return(xx)
		} else {
			if (! inherits(x, 'SpatialPoints') )  {
				x <- data.frame(x[,1:2])
				colnames(x) <- c('x', 'y')
				coordinates(x) <- ~ x + y
			}
			return( .pointsInPolygons(x, object@polygons, sum) )
		}
	}
)

