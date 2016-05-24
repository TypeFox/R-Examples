# Author: Robert J. Hijmans
# Date : December 2009
# Version 1.0
# Licence GPL v3


if (!isGeneric("predict")) {
	setGeneric("predict", function(object, ...)
		standardGeneric("predict"))
}	


.percRank <- function(x, y, tail) {
	b <- apply(y, 1, FUN=function(z)sum(x<z))
	t <- apply(y, 1, FUN=function(z)sum(x==z))
	r <- (b + 0.5 * t)/length(x)
	
	if (tail=='both') {
		i <- which(r > 0.5)
		r[i] <- 1-r[i]
	} else if (tail == 'high') {
		r[ r < 0.5 ] <- 0.5
		r <- 1-r
	} else { # if tail == low
		r[ r > 0.5 ] <- 0.5
	}
	r * 2
}


setMethod('predict', signature(object='Bioclim'), 
function(object, x, tails=NULL, ext=NULL, filename='', useC=TRUE, ...) {

	ln <- colnames(object@presence)

	if (is.null(tails) ) {
		tails <- rep('both', times=length(ln))
	} else {
		test <- all(tails %in% c('low', 'high', 'both'))
		if (!test) {
			stop('"tails" should be a character vector with values "low", "high", "both"')
		}
		if (length(tails) == 1) {
			tails <- rep(tails, times=length(ln))
		} else if (length(tails) != length(ln)) {
			stop('length of "tails" is: ', length(tails), '.\nThis does not match the number of variables in the model which is: ', length(ln))
		}
	}
	tailopt <- match(tails, c('both', 'high', 'low'))
	
	if (! (extends(class(x), 'Raster')) ) {
		if (! all(ln %in% colnames(x)) ) {
			stop('missing variables in x')
		}
		x <- x[, ln ,drop=FALSE]
		
		if (useC) {
			pres <- as.matrix(stats::na.omit(object@presence))
			for (i in 1:ncol(pres)) {
				pres[,i] <- sort(pres[,i])
			}
			if (!inherits(x, 'matrix')) {
				x <- as.matrix(x)
			}
			mincomp <- object@min
			mincomp[tails=='high'] <- -Inf
			maxcomp <- object@max
			maxcomp[tails=='low'] <- Inf
			
			bc <- .Call('percRank', as.double(pres), 
									as.integer(dim(pres)),
									as.double(x),
									as.integer(dim(x)),
									as.double(mincomp), 
									as.double(maxcomp),										 
									as.integer(tailopt), PACKAGE='dismo' )
		} else {
			bc <- matrix(ncol=length(ln), nrow=nrow(x))
			for (i in 1:ncol(bc)) {
				bc[,i] <- .percRank(object@presence[,ln[i]], x[,ln[i], drop=FALSE], tails[i])
			}
			return( apply(bc, 1, min) )
		}
	} else {
		out <- raster(x)
		if (! is.null(ext)) {
			out <- crop(out, ext)
			firstrow <- rowFromY(x, yFromRow(out, 1))
			firstcol <- colFromX(x, xFromCol(out, 1))
		} else {
			firstrow <- 1
			firstcol <- 1
		}
		ncols <- ncol(out)
				
		rasternames <- names(x)
		if (! all(ln %in% rasternames )) {
			stop('missing variables in Raster object')
		}
		if ( inherits(x, 'RasterStack') & (length(ln) < length(rasternames))) {
			x <- x[[ln]]
		}
		
		if (canProcessInMemory(out, 2)) {
			inmem <- TRUE
			v <- matrix(NA, ncol=nrow(out), nrow=ncol(out))
		} else {
			inmem <- FALSE
			out <- writeStart(out, filename, ...)
		}

		tr <- blockSize(out, n=nlayers(x)+2)
		pb <- pbCreate(tr$n, ...)	
		
		mincomp <- object@min
		mincomp[tails=='high'] <- -Inf
		maxcomp <- object@max
		maxcomp[tails=='low'] <- Inf
		
		pres <- as.matrix(stats::na.omit(object@presence))
		for (i in 1:ncol(pres)) {
			pres[,i] <- sort(pres[,i])
		}
		
		for (i in 1:tr$n) {
			rr <- firstrow + tr$row[i] - 1
			vals <- getValuesBlock(x, row=rr, nrows=tr$nrows[i], firstcol, ncols)[, ln, drop=FALSE]
	
			if (useC) {
				res <- .Call('percRank', as.double(pres), 
										 as.integer(dim(pres)),
										 as.double(vals),
										 as.integer(dim(vals)),
										 as.double(mincomp), 
										 as.double(maxcomp),										 
										 as.integer(tailopt), PACKAGE='dismo' )
			
			} else {
				bc <- matrix(0, ncol=ncol(vals), nrow=nrow(vals))
				na <- as.vector( attr(stats::na.omit(vals), 'na.action') )
				bc[na] <- NA
				k <- (apply(t(vals) >= mincomp, 2, all) & apply(t(vals) <= maxcomp, 2, all))
				k[is.na(k)] <- FALSE
				for (j in 1:length(ln)) {
					bc[k,j] <- .percRank( pres[ ,ln[j]], vals[k, ln[j], drop=FALSE], tails[j] )
				}
				res <- apply(bc, 1, min)
			}
			
			
			if (inmem) {
				res <- matrix(res, nrow=ncols)
				cols <- tr$row[i]:(tr$row[i]+dim(res)[2]-1)
				v[ , cols] <- res
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
}
)


