# Author: Robert J. Hijmans
# contact: r.hijmans@gmail.com
# Date : December 2009
# Version 0.1
# Licence GPL v3

setMethod('predict', signature(object='Mahalanobis'), 
function(object, x, toCenter=FALSE, ext=NULL, filename='', ...) {

	if (toCenter) {
		m <- colMeans(object@presence, na.rm=TRUE)
	}

	if (! (extends(class(x), 'Raster')) ) {
		if (! all(colnames(object@presence) %in% colnames(x)) ) {
			stop('missing variables in matrix ')
		}
		x <- x[ , colnames(object@presence),drop=FALSE]
		if (toCenter) {
			mah <- 1 - apply(data.frame(x), 1, FUN=function(z) mahalanobis(z, m, object@cov, inverted=TRUE))
		return(mah)
		} else {
			mah <- 1 - apply(data.frame(x), 1, FUN=function(z) min( mahalanobis(object@presence, z, object@cov, inverted=TRUE)))
		}
	} else {
	
		if (! all(colnames(object@presence) %in% names(x)) ) {
			stop('missing variables in Raster object ')
		}
	
		out <- raster(x)
		if (!is.null(ext)) {
			out <- crop(out, ext)
			firstrow <- rowFromY(x, yFromRow(out, 1))
			firstcol <- colFromX(x, xFromCol(out, 1))
		} else {
			firstrow <- 1
			firstcol <- 1
		}
		ncols <- ncol(out)

			
		if (canProcessInMemory(out, 2)) {
			inmem <- TRUE
			v <- matrix(NA, ncol=nrow(out), nrow=ncol(out))
		} else {
			inmem <- FALSE
			if  (filename == '') {
				filename <- rasterTmpFile()
				if (getOption('verbose')) { message('writing raster to:', filename)	}						
			}
			out <- writeStart(out, filename=filename, ...)
		}

		cn <- colnames(object@presence)
		
		tr <- blockSize(out, n=nlayers(x)+2)
		pb <- pbCreate(tr$n, ...)	
		for (i in 1:tr$n) {
			rr <- firstrow + tr$row[i] - 1
			vals <- getValuesBlock(x, row=rr, nrows=tr$nrows[i], firstcol, ncols)

			vals <- vals[,cn,drop=FALSE]
			if (toCenter) {
				res <- 1 - apply(data.frame(vals), 1, FUN=function(z) mahalanobis(z, m, object@cov, inverted=TRUE))
			} else {
				res <- 1 - apply(data.frame(vals), 1, FUN=function(z) min( mahalanobis(object@presence, z, object@cov, inverted=TRUE)))
			}
			if (inmem) {
				res <- matrix(res, nrow=ncol(out))
				cols = tr$row[i]:(tr$row[i]+dim(res)[2]-1)
				v[, cols] <- res
			} else {
				out <- writeValues(out, res, tr$row[i])
			}
			pbStep(pb, i) 

		} 
		if (inmem) {
			out <- setValues(out, as.vector(v))
			if (filename != '') {
				out <- writeRaster(out, filename, ...)
			}
		} else {
			out <- writeStop(out)	
		}

		pbClose(pb)
		return(out)
	}
})


