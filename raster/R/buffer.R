# Author: Robert J. Hijmans
# Date : September 2009
# Version 0.9
# Licence GPL v3


if (!isGeneric('buffer')) {
	setGeneric('buffer', function(x, ...)
		standardGeneric('buffer'))
}	

setMethod('buffer', signature(x='Spatial'), 
function(x, width=1, dissolve=TRUE, ...) {
	stopifnot(requireNamespace("rgeos"))
	rgeos::gBuffer(x, byid=!dissolve, width=width, ...)
}
)


setMethod('buffer', signature(x='RasterLayer'), 
function(x, width=0, filename='', doEdge=FALSE, ...) {

	stopifnot(width > 0)

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
		pb <- pbCreate(4, label='buffer', ...)
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
		x[x > width] <- NA
		x[!is.na(x)] <- 1
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
	pb <- pbCreate(tr$n, label='buffer', ...)
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
		vals[vals > width] <- NA
		vals[!is.na(vals)] <- 1
		out <- writeValues(out, vals, tr$row[i])
		pbStep(pb) 	
	}	
	pbClose(pb)
	out <- writeStop(out)
	return(out)
}
)

