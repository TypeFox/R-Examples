# Author: Robert J. Hijmans
# Date : September 2009
# Version 0.9
# Licence GPL v3

if (!isGeneric("distance")) {
	setGeneric("distance", function(x, ...)
		standardGeneric("distance"))
}	


setMethod('distance', signature(x='RasterLayer'), 
function(x, filename='', doEdge=TRUE, ...) {
	if (doEdge) {
		r <- boundaries(x, classes=FALSE, type='inner', progress=.progress(...)) 
		pts <- try(  rasterToPoints(r, fun=function(z){ z>0 } )[,1:2, drop=FALSE] )
	} else {
		pts <- try(  rasterToPoints(x)[,1:2, drop=FALSE] )
	}
	
	if (class(pts) == "try-error") {
		return( .distanceRows(x, filename=filename, ...) )
	}
	if (nrow(pts) == 0) {
		stop('RasterLayer has no NA cells (for which to compute a distance)')
	}
	out <- raster(x)
	filename <- trim(filename)
	
	if (couldBeLonLat(x)) { 
		longlat=TRUE 
	} else { 
		longlat=FALSE 
	}
	                                                                        
	if (canProcessInMemory(out, 6)) {
		pb <- pbCreate(3, label='distance', ...)
		x <- values(x)
		i <- which(is.na(x))
		if (length(i) < 1) {
			stop('raster has no NA values to compute distance to')
		}
		pbStep(pb)
		x[] <- 0
		xy <- xyFromCell(out, i)
		x[i] <- .Call("distanceToNearestPoint", xy, pts, as.integer(longlat), PACKAGE='raster')
		pbStep(pb)
		out <- setValues(out, x)
		if (filename != '') {
			out <- writeRaster(out, filename=filename, ...)
		}
		pbStep(pb)
		pbClose(pb)
		return(out)
	} 
	
	out <- writeStart(out, filename=filename, ...)
	tr <- blockSize(out)
	pb <- pbCreate(tr$n, label='distance', ...)
	xy <- cbind(rep(xFromCol(out, 1:ncol(out)), tr$nrows[1]), NA)
	for (i in 1:tr$n) {
		if (i == tr$n) {
			xy <- xy[1:(ncol(out)*tr$nrows[i]), ]
		}
		xy[,2] <- rep(yFromRow(out, tr$row[i]:(tr$row[i]+tr$nrows[i]-1)), each=ncol(out))
		vals <- getValues(x, tr$row[i], tr$nrows[i])
		j <- which(is.na(vals))
		vals[] <- 0
		if (length(j) > 0) {
			vals[j] <- .Call("distanceToNearestPoint", xy[j,,drop=FALSE], pts, as.integer(longlat), PACKAGE='raster')
		}
		out <- writeValues(out, vals, tr$row[i])
		pbStep(pb) 	
	}	
	pbClose(pb)
	out <- writeStop(out)
	return(out)
}
)

