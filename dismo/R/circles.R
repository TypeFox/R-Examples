# Author: Robert J. Hijmans
# Date : Febrary 2010
# Version 0.1
# Licence GPL v3


.avgDist <- function(xy, lonlat, r=6378137) {
	xy <- unique(xy)
	xy <- as.matrix(xy[,1:2])
	d <- matrix(nrow=nrow(xy), ncol=nrow(xy))
	if (lonlat) {
		for (i in 1:nrow(xy)) { 
			d[,i] <- pointDistance(xy[i,], xy, longlat=TRUE, r=r)
		}
	} else {
		for (i in 1:nrow(xy)) { 
			d[,i] <- pointDistance(xy[i,], xy, longlat=FALSE, r=r)
		}
	}
	median(apply(d, 1, function(x)quantile(x, 0.1)))
}


.generateCircles <- function(xy, d, n=360, lonlat, r=6378137, crs=NA ) {
	if (length(d)==1) {
		d <- rep(d, nrow(xy))
	} else if (length(d) != nrow(xy)) {
		# recycling
		dd <- vector(length=nrow(xy))
		dd[] <- d
		d <- dd
	}
	
	xy <- as.matrix(xy[,1:2])
	n <- max(4, round(n))
	toRad <- pi/180
	brng <- 1:n * 360/n
	brng <- brng * toRad
	if (lonlat) { 
		xy <- xy * toRad 
	}
	pols <- list()
	for (i in 1:nrow(xy)) {
		p <- xy[i, ]
		if (lonlat) {
			lon1 <- p[1]
			lat1 <- p[2]
			lat2 <- asin(sin(lat1) * cos(d[i]/r) + cos(lat1) * sin(d[i]/r) * cos(brng))
			lon2 <- lon1 + atan2(sin(brng) * sin(d[i]/r) * cos(lat1), cos(d[i]/r) - sin(lat1) * sin(lat2))
			lon2 <- (lon2 + pi)%%(2 * pi) - pi
			lon2[is.nan(lon2)] <- NA
			lat2[is.nan(lat2)] <- NA
			res <- cbind(lon2, lat2)/toRad
			colnames(res) <- c("lon", "lat")
		} else {
			x2 <- p[1] + d[i] * cos(brng)
			y2 <- p[2] + d[i] * sin(brng)
			res <- cbind(x2, y2)
			colnames(res) <- c("x", "y")
		}
		res <- rbind(res, res[1,])
		pols <- c(pols, Polygons(list(Polygon( res )), i))		
	}

	pols <- SpatialPolygons(pols)
	projection(pols) <- crs
	return( pols )
}


setClass('CirclesRange',
	contains = 'DistModel',
	representation (
		polygons='SpatialPolygons'
	),	
	prototype (	
	),
	validity = function(object)	{
		return(TRUE)
	}
)


setMethod("geometry", "CirclesRange",
	function(obj) {
		obj@polygons
	}
)

setMethod("polygons", "CirclesRange",
	function(obj) {
		obj@polygons
	}
)


if (!isGeneric("circles")) {
	setGeneric("circles", function(p, ...)
		standardGeneric("circles"))
}	

setMethod('circles', signature(p='data.frame'), 
	function(p, d, lonlat, n=360, r=6378137, dissolve=TRUE, ...) {
		ci <- new('CirclesRange')
		ci@presence <- p
		if (missing(lonlat)) {
			stop("you must provide a 'lonlat' argument")
		}
		if (missing(d)) {
			d <- .avgDist(p, lonlat=lonlat, ...) / 2
		}
		
		if (inherits(p, 'SpatialPoints')) {
			crs <- projection(p)
		} else if (lonlat) {
			crs <- '+proj=longlat +datum=WGS84'		
		} else {
			crs <- NA
		}
		
		ci@polygons <- .generateCircles(p, d=d, lonlat=lonlat, crs=crs, n=n, r=r, ...)
		if (dissolve & requireNamespace('rgeos')) {
			ci@polygons <- rgeos::gUnaryUnion(ci@polygons)
		}
		return(ci)
	}
)


setMethod('circles', signature(p='matrix'), 
	function(p, d, lonlat, n=360, r=6378137, dissolve=TRUE, ...) {
		if (missing(d)) { 
			d <- .avgDist(p, lonlat=lonlat, r=r) / 2 
		}
		circles(data.frame(p), d=d, lonlat=lonlat, n=n, r=r, ...)
	}
)

setMethod('circles', signature(p='SpatialPoints'), 
	function(p, d, lonlat, n=360, r=6378137, dissolve=TRUE, ...) {
		if (missing(lonlat)) {
			lonlat <- isLonLat(p)
		}
		p <- coordinates(p)
		if (missing(d)) { 
			d <- .avgDist(p, lonlat=lonlat, r=r) / 2	
		}
		circles(p, d=d, lonlat=lonlat, n=n, r=r, ...)
	}
)


setMethod("plot", signature(x='CirclesRange', y='missing'), 
	function(x, ...) {
		sp::plot(x@polygons, ...)
	}
)

setAs('CirclesRange', 'SpatialPolygons', 
	function(from) {
		from@polygons
	}
)


