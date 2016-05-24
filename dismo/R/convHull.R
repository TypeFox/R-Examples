# Author: Robert J. Hijmans
# Date : Febrary 2010
# Version 0.1
# Licence GPL v3


setClass('ConvexHull',
	contains = 'DistModel',
	representation (
		polygons='SpatialPolygonsDataFrame'
	),	
	prototype (	
	),
	validity = function(object)	{
		return(TRUE)
	}
)


setMethod("polygons", "ConvexHull",
	function(obj) {
		obj@polygons
	}
)

if (!isGeneric("convHull")) {
	setGeneric("convHull", function(p, ...)
		standardGeneric("convHull"))
}	


setMethod('convHull', signature(p='data.frame'), 
	function(p, n=1, crs=NULL, ...) {
		ch <- new('ConvexHull')
		ch@presence <- p
		ch@polygons <- .generateHulls(p, n)
		if (!is.null(crs)) {
			projection(ch@polygons) <- crs
		}
		return(ch)
	}
)


setMethod('convHull', signature(p='matrix'), 
	function(p, ...) {
		convHull(as.data.frame(p), ...)
	}
)

setMethod('convHull', signature(p='SpatialPoints'), 
	function(p, ...) {
		convHull(coordinates(p), crs=p@proj4string, ...)
	}
)



.generateHulls <- function(xy, n=1) {
	xy <- unique(  stats::na.omit(xy[, 1:2]) )
    if (nrow(xy) < 3) { stop ('Insuficient number of points to make a hull; you need at least 3 unique points' ) }
    n <- pmax(1, round(n))
    n <- pmin(n, floor(nrow(xy) / 3))
    n = unique(n)
    pols = list()

    count <- 1
    for (k in n) {
		if (k == 1) {
			h <- xy[chull(xy), ]
			pols <- c(pols, Polygons(list(Polygon( rbind(h, h[1,]) )), count))
		} else {
			ch = integer()
			cl = kmeans(xy, k, 100)$cluster
			clusters = unique(cl)
			subp = list()
			for (i in clusters) {
				pts <- xy[cl==i, ]
				h <- pts[chull(pts), ]
				subp <- c(subp, Polygon( rbind(h, h[1,]) ))
			}
			pols <- c(pols, Polygons( subp, count) )
		}
		count <- count + 1
	}
	pols <- SpatialPolygonsDataFrame(SpatialPolygons( pols ), data.frame(id=n, w=1/length(n)) )
    return( pols )
}


