# Author: Robert J. Hijmans
# Date : Febrary 2010
# Version 0.1
# Licence GPL v3


.pointsInPolygons <- function(xy, polygons, fun=NULL) {
# based on similar function in sp (overlay) by Pebesma and Bivand
# and calling a C function from the SP package
# the difference is that this functions considers that a point can be inside
# multiple polygons
	polygons <- polygons@polygons
	xy <- coordinates(xy)
	xy[] <- as.numeric(xy)
	res <- vector(length=nrow(xy))
	result <- matrix(nrow=length(res), ncol=length(polygons))
	for (i in 1:length(polygons)) {
		res[] <- 0
		poly <- polygons[[i]]
		e <- extent(bbox(poly))
		inbox <- xy[,1] >= e@xmin & xy[,1] <= e@xmax & xy[,2] >= e@ymin & xy[,2] <= e@ymax
		p <- xy[inbox, ,drop=FALSE]
		if (nrow(p) > 0) {
			res2 <- vector(length=nrow(p))
			res2[] <- 0
			for (j in 1:length(poly@Polygons)) {
				pol.x <- poly@Polygons[[j]]@coords[,1]
				pol.y <- poly@Polygons[[j]]@coords[,2]
				resj <- point.in.polygon(p[,1], p[,2], pol.x, pol.y)
				#resj <- .Call("R_point_in_polygon_sp", as.numeric(p[,1]), as.numeric(p[,2]), as.numeric(pol.x), as.numeric(pol.y), PACKAGE = "sp")
				resj <- as.logical(resj)
				if (poly@Polygons[[j]]@hole) {
					resj[resj == 1] <- NA
				}
				res2 <- pmax(res2, resj)
			}
			res[inbox] <- pmax(res[inbox], res2)	
		}
		result[,i] <- res
	}
	result[is.na(result)] = 0
	if (! is.null(fun)) {
		result <- apply(result, 1, fun)
	}
	return(result)
}


# a =  pointInPolygon(xy, wrld_simpl)
# b = apply(a, 1, max)
# c = apply(a, 1, which.max)
# c[b==0] = NA


#..xpointsInPolygons <- function(xy, polygons, fun=NULL) {
#	stopifnot( require(rgeos) )
	#proj4string(xy) <- proj4string(polygons)
#	result <- gIntersects(xy, polygons, byid=TRUE)
#	if (! is.null(fun)) {
#		result <- apply(result, 1, fun)
#	}
#	result * 1
#}


# a =  pointInPolygon(xy, wrld_simpl)
# b = apply(a, 1, max)
# c = apply(a, 1, which.max)
# c[b==0] = NA

