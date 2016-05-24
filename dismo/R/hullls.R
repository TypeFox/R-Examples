# Author: Robert J. Hijmans
# contact: r.hijmans@gmail.com
# Date : Febrary 2010
# Version 0.1
# Licence GPL v3


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


#hullmodel <- function(xy) {
#       n = max(1, log(nrow(xy)))
#       return(  hulls(xy, 1:n) )
#}

#predict.hull <- function(r, hullmodel, ...) {
#	rasterize(hullmodel, r, field=1, overlap='sum', ...)
#}
