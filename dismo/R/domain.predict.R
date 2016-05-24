# Author: Robert J. Hijmans
# Date : December 2009
# Version 0.1
# Licence GPL v3


setMethod('predict', signature(object='Domain'), 
function(object, x, ext=NULL, filename='', ...) {

	domdist <- function(xx, ii, y) {
		r <- xx@range[ii]
		xx <- xx@presence[,ii]
		# improvement by David Woolford
		d <- apply(data.frame(y), 1, FUN=function(z)(mean(abs(xx-z)/r)))
		#d <- apply(data.frame(y), 1, FUN=function(z)(abs(xx-z)/r))
		#d <- apply(d, 2, mean)
		d[which(d > 1)] <- 1
		1-d
	}
	
	domdistcat <- function(xx, ii, y) {
		stop("can't do factors yet")
	}
		
	if (! (extends(class(x), 'Raster')) ) {
		if (! all(colnames(object@presence) %in% colnames(x)) ) {
			stop('missing variables in x ')
		}
		
		dom <- matrix(ncol=length(colnames(object@presence)), nrow=nrow(x) )
		ln <- colnames(object@presence)
		
		f <- which(ln %in% object@factors)
		if (length(f) > 0 ) {
			i1 <- 1:ncol(dom)[-f]
			i2 <- 1:ncol(dom)[f]
			for (i in i1) {
				dom[,i] <- domdist(object, ln[i], x[,ln[i]])
			}
			for (i in i2) {
				dom[,i] <- domdistcat(object, ln[i], x[,ln[i]])
			}
		} else {
			for (i in 1:ncol(dom)) {
				dom[,i] <- domdist(object, ln[i], x[,ln[i]])
			}
		}
		return ( apply(dom, 1, min ) )

	} else {

		if (! all(colnames(object@presence) %in% names(x)) ) {
			stop('missing variables in Raster object')
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
			inmem=TRUE
			v <- matrix(NA, ncol=nrow(out), nrow=ncol(out))
		} else {
			inmem <- FALSE
			if  (filename == '') {
				filename <- rasterTmpFile()
				if (getOption('verbose')) { message('writing raster to:', filename)	}						
			}
			out <- writeStart(out, filename=filename, ... )
		}
		
		ln <- colnames(object@presence)
		tr <- blockSize(out, n=nlayers(x)+2)
		dom <- matrix(ncol=nlayers(x), nrow=ncols*tr$nrows[1] )
		pb <- pbCreate(tr$n, ...)	
		
		for (i in 1:tr$n) {
			if (i == tr$n) {
				dom <- matrix(ncol=nlayers(x), nrow=ncols*tr$nrows[i] )
			}
		
			rr <- firstrow + tr$row[i] - 1
			vals <- getValuesBlock(x, row=rr, nrows=tr$nrows[i], firstcol, ncols)

			for (j in 1:ncol(dom)) {
				dom[,j] <- domdist(object, ln[j], vals[,ln[j]])
			}
			res <- apply(dom, 1, min)
			if (inmem) {
				v[, tr$row[i]:(tr$row[i]+tr$nrows[i]-1)] <- matrix(res, nrow=ncol(out))
			} else {
				out <- writeValues(out, res, tr$row[i])
			}
			pbStep(pb, i) 
		} 
		
		if (inmem) {
			out <- setValues(out, as.vector(v))
			if (filename != '') {
				out <- writeRaster(out, filename=filename, ...)
			}
		} else {
			out <- writeStop(out)	
		}
		pbClose(pb)
		return(out)
	}
}

)

