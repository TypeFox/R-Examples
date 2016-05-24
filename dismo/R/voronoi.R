# Author: Robert J. Hijmans
# Date : Febrary 2010
# Version 1.0
# Licence GPL v3

# adapted from code by Carson Farmer
# http://www.carsonfarmer.com/?p=455
voronoi <- function(xy){

	if (!requireNamespace('deldir')) { stop('you need to first install the deldir libary') }

	dat <- NULL
	sp <- FALSE
	if (inherits(xy, 'Spatial')) {
		if (.hasSlot(xy, 'data')) {
			dat <- slot(xy, 'data')
		}
		prj <- proj4string(xy)
		sp <- TRUE
		xy <- coordinates(xy)
		dups <- duplicated(xy)
		if (any(dups)) {
			xy <- xy[!dups, ,drop=FALSE]
			dat <- dat[!dups, ,drop=FALSE]
		}
	} else {
		xy <- stats::na.omit(xy[, 1:2])
		xy <- unique(xy)
	}
	
	z <- deldir::deldir(xy[,1], xy[,2])
	w <- deldir::tile.list(z)
	polys <- vector(mode='list', length=length(w))
	for (i in seq(along=polys)) {
		pcrds <- cbind(w[[i]]$x, w[[i]]$y)
		pcrds <- rbind(pcrds, pcrds[1,])
		polys[[i]] <- Polygons(list(Polygon(pcrds)), as.character(i))
	}
	if (sp) {
		polys <- SpatialPolygons(polys, proj4string=CRS(prj))
	} else {
		polys <- SpatialPolygons(polys)
	}
	
	if (is.null(dat)) {
		dat <- data.frame(id=1:nrow(xy), xy)
	}
	rownames(dat) <- row.names(polys)
	polys <- SpatialPolygonsDataFrame(polys, data=dat)
	return(polys)
}

