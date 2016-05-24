# Author: Robert J. Hijmans
# Date : November 2012
# Version 1.0
# Licence GPL v3


setClass('EcoLim',
	contains = 'DistModel',
	representation (
		funs = 'list',
		x = 'matrix',
		y = 'matrix'
	),	
	prototype (	
	),
	validity = function(object)	{
		return(TRUE)
	}
)


if (!isGeneric("ecolim")) {
	setGeneric("ecolim", function(x, y, ...)
		standardGeneric("ecolim"))
}	

setMethod('ecolim', signature(x='matrix', y='matrix'), 
	function(x, y=matrix(c(0,0,1,1,0,0)), extrapolate=TRUE,...) {
		stopifnot(NCOL(x) == NCOL(y))
		cn <- colnames(x)
		if (any(cn == "")) {
			stop("All columns of 'x' must have names")
		}
		m <- new('EcoLim')
		f <- list()
		extrapolate <- as.logical(extrapolate)+1
		extrapolate <- rep(extrapolate, length.out=ncol(x))
		for (i in 1:ncol(x)) {
			xy <- stats::na.omit(cbind(x[,i], y[,i]))
			f[[i]] <- approxfun(xy[,1], xy[,2], rule=extrapolate[i],...)
		}
		m@funs <- f
		m@x <- x
		m@y <- y
		m
	}
)


setMethod('show', signature(object='EcoLim'), 
function(object) {
	cat('class     : EcoLim\n')
	cat('variables :', colnames(object@x), '\n')
} )


setMethod('plot', signature(x='EcoLim'), 
function(x, ...) {
	n <- ncol(x@x)
	nc <- floor(sqrt(n))
	nr <- ceiling(sqrt(n))
	graphics::par(mfrow=c(nr,nc))
	nm <- colnames(x@x)
	for (i in 1:n) {
		plot(x@x[,i], x@funs[[i]](x@x[,i]), type='l', xlab=nm[i], ylab='response', ...)
	}
} )


setMethod('predict', signature(object='EcoLim'), 
function(object, x, fun=min, ext=NULL, filename='', ...) {

	fn <- colnames(object@x)
	
	if (! (extends(class(x), 'Raster')) ) {
		ln <- colnames(x)
		if (! all(fn %in% ln) ) {
			stop('missing variables in x')
		}
		x <- x[, fn ,drop=FALSE]
		ec <- matrix(ncol=length(fn), nrow=nrow(x))
		for (i in 1:ncol(ec)) {
			ec[,i] <- object@funs[[i]](x[,i])
		}
		return( apply(ec, 1, fun) )
		
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
				
		ln <- names(x)
		if (! all(fn %in% ln )) {
			stop('missing variables in Raster object')
		}
		if ( (inherits(x, 'RasterStack')) & (length(fn) < length(ln))) {
			x <- x[[fn]]
		}
		
		if (canProcessInMemory(out, 2)) {
			inmem=TRUE
			v <- matrix(NA, ncol=nrow(out), nrow=ncol(out))
		} else {
			inmem <- FALSE
			if  (filename == '') {
				filename <- rasterTmpFile()
				out <- writeStart(out, filename, ...)
			}
		}

		tr <- blockSize(out, n=nlayers(x)+2)
		pb <- pbCreate(tr$n, ...)	
		
		for (i in 1:tr$n) {
			rr <- firstrow + tr$row[i] - 1
			vals <- getValuesBlock(x, row=rr, nrows=tr$nrows[i], firstcol, ncols)[, fn, drop=FALSE]
	
			ec <- matrix(0, ncol=ncol(vals), nrow=nrow(vals))
			for (j in 1:length(fn)) {
				ec[,j] <- object@funs[[j]](vals[,j])
			}
			res <- apply(ec, 1, fun)
			
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

