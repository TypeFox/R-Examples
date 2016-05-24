# Author: Robert J. Hijmans
# contact: r.hijmans@gmail.com
# Date : Febrary 2010
# Version 0.1
# Licence GPL v3


setMethod('predict', signature(object='GeographicDistance'), 
	function(object, x, ext=NULL, mask=FALSE, scale=1, fun=NULL, filename='', ...) {
			
		if (is.null(fun)) {
			inverse <- function(x, ...) {
				x <- x / scale
				x[x < 1] <- 1
				1/x
			}
		} else {
			inverse <- function(x, ...) {
				x <- x / scale
				x[x < 1] <- 1
				fun(1/x)
			}
		}

		if ( extends(class(x), 'Raster'))  {
			if (!hasValues(x)) {
				mask <- FALSE
			}
			if (! mask) {
				x <- raster(x)
			}
			if (! is.null(ext)) { 
				x <- crop(x, ext) 
			}

			xx <- distanceFromPoints(x, object@presence)
			if (mask) {
				xx <- mask(xx, x)
			}
			xx <- calc(xx, fun=inverse, filename=filename, ...)
			return(xx)
			
		} else {
		
			if ( inherits(x, 'SpatialPoints') )  { 
				x <- coordinates(x) 
			} else {
				x <- x[, colnames(object@presence)]
			}
			
			res <- vector(length=nrow(x))
			for (i in 1:nrow(x)) {
				res[i] <- min( pointDistance(x[i,], object@presence, longlat=object@lonlat) )
			}	
			return(inverse(res))
		}
	}
)

