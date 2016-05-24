# Author: Robert J. Hijmans
# contact: r.hijmans@gmail.com
# Date : April 2010
# Version 0.1
# Licence GPL v3


setMethod('predict', signature(object='VoronoiHull'), 
	function(object, x, ext=NULL, filename='', mask=FALSE, ...) {
	
		if ( inherits(x, 'Raster'))  {
			if (! mask) {
				x <- raster(x)
			}
			if (! is.null(ext)) { 
				x <- crop(x, ext) 
			}
			
			if (mask) {
				xx <- rasterize(object@polygons, raster(x), field='pa', mask=FALSE, update=FALSE, getCover=FALSE, silent=TRUE, ...)
				xx <- mask(xx, x, filename=filename, ...)
			} else {
				xx <- rasterize(object@polygons, raster(x), field='pa', mask=FALSE, update=FALSE, getCover=FALSE, silent=TRUE, filename=filename, ...)
			}
			
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

