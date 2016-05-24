# Author: Robert J. Hijmans
# contact: r.hijmans@gmail.com
# Date : Febrary 2010
# Version 0.1
# Licence GPL v3


setClass('GeographicDistance',
	contains = 'DistModel',
	representation (
		lonlat='logical'
	),	
	prototype (	
	),
	validity = function(object)	{
		return(TRUE)
	}
)


if (!isGeneric("geoDist")) {
	setGeneric("geoDist", function(p, ...)
		standardGeneric("geoDist"))
}	


setMethod('geoDist', signature(p='data.frame'), 
	function(p, lonlat, ...) {
		gd <- new('GeographicDistance')
		gd@presence <- stats::na.omit(p)
		gd@lonlat <- lonlat
		return(gd)
	}
)


setMethod('geoDist', signature(p='matrix'), 
	function(p, lonlat, ...) {
		stopifnot(ncol(p) == 2)
		p <- as.data.frame(p)
		geoDist(p, lonlat=lonlat, ...) 
	}
)

setMethod('geoDist', signature(p='SpatialPoints'), 
	function(p, lonlat, ...) {
		if (missing(lonlat)) {
			lonlat <- isLonLat(p)
		}
		p <- coordinates(p)
		geoDist(p, lonlat=lonlat, ...) 
	}
)

