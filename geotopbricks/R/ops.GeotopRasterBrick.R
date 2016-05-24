


NULL
#' 
#' \code{Ops} method for a \code{GeotopRasterBrick} object
#' 
#' @param e1,e2  the \code{\link{GeotopRasterBrick}} or numeric  objects
#' @note If \code{e1} or \code{e2} time index is not taken into account. 
#' 
#' @title Ops
#' @name Ops
#' @importFrom "raster" subset
#' @importFrom "methods" Arith 
#' @importFrom "methods" Compare
#' @importFrom "methods" Logic
#' @export
#' @rdname Ops-methods
#' @keywords methods
#' @docType methods
#' @method Ops GeotopRasterBrick GeotopRasterBrick
#' @aliases Ops,GeotopRasterBrick,GeotopRasterBrick-method



setMethod("Ops", signature(e1='GeotopRasterBrick', e2='GeotopRasterBrick'),
		function(e1,e2) {
		
		
		
		out <- new('GeotopRasterBrick')	
			
	
		index <- e1@index[e1@index %in% e2@index]
		
		
		i1 <- which(sort(e1@index) %in% index)
		i2 <- which(sort(e2@index) %in% index)
		
		b1a <- brick(e1)
		b2a <- brick(e2)
#		if (is(b1a,"RasterBrick")){
#			
#			return(b1a)
#		}
		
		b1 <- subset(b1a,i1)
		b2 <- subset(b2a,i2)
		
		out@ascpath <- as.zoo(array(NA,length(index)))
		index(out@ascpath) <- index 
		
		
		
		out@index <- index
		
		
		
		out@brick <- callGeneric(b1,b2)
# These commnted code lines do not work!!! 		
#		if ((e1@layer)==(e2@layer)) {
#			out@layer <- e1@layer 
#		} else {
#			
#			out@layer <- as.character(NA) 
#		}
#	END !! THE BUG HAS NOT YET BEEN SOLVED!! 		
		out@layer <- e1@layer
		
		return(out)
			
		}
)

NULL

#'
#' @export
#' @rdname Ops-methods
#' 
#' 
#' @method Ops GeotopRasterBrick numeric
#' @aliases Ops,GeotopRasterBrick,numeric-method
setMethod("Ops", signature(e1='GeotopRasterBrick', e2='numeric'),
		function(e1,e2) {
			
			out <- e1
			
			
			out@brick <- callGeneric(brick(e1),e2)
			
			return(out)
			
		}
)

NULL

#'
#' @export
#' @rdname Ops-methods
#' 
#' 
#' 
#' @rdname Ops-methods
# @keywords methods
# @docType methods
#' @method Ops numeric GeotopRasterBrick
#' @aliases Ops,numeric,GeotopRasterBrick-method

setMethod("Ops", signature(e1='numeric', e2='GeotopRasterBrick'),
		function(e1,e2) {
			
			out <- e2
			
			
			out@brick <- callGeneric(e1,brick(e2))
			
			return(out)
			
		}
)

#
#setMethod("Arith", signature(e1='Raster', e2='Raster'),
#    function(e1, e2){ 
#
#		if (!hasValues(e1)) { stop('first Raster object has no values') }
#		if (!hasValues(e2)) { stop('second Raster object has no values') }
#		
#		nl1 <- nlayers(e1)
#		nl2 <- nlayers(e2)
#		nl <- max(nl1, nl2)
#
#		proj1 <- projection(e1)
#		proj2 <- projection(e2)
#	
#		if ( ! compare(e1, e2, crs=FALSE, stopiffalse=FALSE) ) {
#			if ( compare(e1, e2, extent=FALSE, rowcol=FALSE, crs=TRUE, res=TRUE, orig=TRUE, stopiffalse=TRUE) ) {
#				ie <- intersect(extent(e1), extent(e2))
#				if (is.null(ie)) { 	stop() }
#				warning('Raster objects have different extents. Result for their intersection is returned')
#				e1 <- crop(e1, ie)
#				e2 <- crop(e2, ie)
#			} else {
#				stop()  # stops anyway because compare returned FALSE
#			}
#		}
#
#		if (nl > 1) {
#			r <- brick(e1, values=FALSE, nl=nl)
#		} else {
#			r <- raster(e1)
#		}
#
#		
#		if (canProcessInMemory(r, 4)) {
#			if (nl1 == nl2 ) {
#				return( setValues(r, values=callGeneric( getValues(e1), getValues(e2))) )
#			} else {
#				return( setValues(r, matrix(callGeneric( as.vector(getValues(e1)), as.vector(getValues(e2))), ncol=nl)) )
#			}
#			
#		} else {
#		
#			tr <- blockSize(e1)
#			pb <- pbCreate(tr$n)			
#			r <- writeStart(r, filename=rasterTmpFile(), overwrite=TRUE )
#			if (nl1 == nl2 ) {
#				for (i in 1:tr$n) {
#					v1 <- getValues(e1, row=tr$row[i], nrows=tr$nrows[i])
#					v2 <- getValues(e2, row=tr$row[i], nrows=tr$nrows[i])
#					v <- callGeneric( v1, v2 )
#					r <- writeValues(r, v, tr$row[i])
#					pbStep(pb, i) 	
#				}
#			} else {
#				for (i in 1:tr$n) {
#					v1 <- as.vector(getValues(e1, row=tr$row[i], nrows=tr$nrows[i]))
#					v2 <- as.vector(getValues(e2, row=tr$row[i], nrows=tr$nrows[i]))
#					v <- matrix(callGeneric( v1, v2 ), ncol=nl)
#					r <- writeValues(r, v, tr$row[i])
#					pbStep(pb, i) 	
#				}
#			}
#			r <- writeStop(r)
#			pbClose(pb)
#			return(r)
#			
#		}
#	}	
#)
#
#
#
#
#setMethod("Arith", signature(e1='RasterLayer', e2='numeric'),
#    function(e1, e2){ 
#		if (!hasValues(e1)) { stop('RasterLayer has no values') }
#
#		r <- raster(e1)
#		if (canProcessInMemory(e1, 4)) {
#			if (length(e2) > ncell(r)) {
#				e2 <- e2[1:ncell(r)]
#			}
#			return ( setValues(r,  callGeneric(as.numeric(getValues(e1)), e2) ) )
#			
#		} else {
#			tr <- blockSize(e1)
#			pb <- pbCreate(tr$n)			
#			r <- writeStart(r, filename=rasterTmpFile(), format=.filetype(), overwrite=TRUE )
#
#			if (length(e2) > 0) {
#				for (i in 1:tr$n) {
#					e <- .getAdjustedE(r, tr, i, e2)
#					v <- callGeneric(getValues(e1, row=tr$row[i], nrows=tr$nrows[i]), e)
#					r <- writeValues(r, v, tr$row[i])
#					pbStep(pb, i) 	
#				}
#			} else {
#				for (i in 1:tr$n) {
#					v <- callGeneric( getValues(e1, row=tr$row[i], nrows=tr$nrows[i]), e2 )
#					r <- writeValues(r, v, tr$row[i])
#					pbStep(pb, i)
#				}
#			}
#			r <- writeStop(r)
#			pbClose(pb)
#			return(r)
#		}		
#	}
#)
#
#
#
#setMethod("Arith", signature(e1='numeric', e2='RasterLayer'),
#    function(e1, e2){ 
#		stopifnot(hasValues(e2))
#
#		r <- raster(e2)
#		if (canProcessInMemory(e2, 4)) {
#			if (length(e1) > ncell(r)) {
#				e1 <- e1[1:ncell(r)]
#			}
#			return ( setValues(r,  callGeneric(e1, getValues(e2)) ) )
#			
#		} else {
#			tr <- blockSize(e2)
#			pb <- pbCreate(tr$n)			
#			r <- writeStart(r, filename=rasterTmpFile(), format=.filetype(), overwrite=TRUE )
#
#			if (length(e1) > 0) {
#				for (i in 1:tr$n) {
#					e <- .getAdjustedE(r, tr, i, e1)
#					v <- callGeneric(e, getValues(e2, row=tr$row[i], nrows=tr$nrows[i]))
#					r <- writeValues(r, v, tr$row[i])
#					pbStep(pb, i) 	
#				}
#			} else {
#				for (i in 1:tr$n) {
#					v <- callGeneric(e1, getValues(e2, row=tr$row[i], nrows=tr$nrows[i]))
#					r <- writeValues(r, v, tr$row[i])
#					pbStep(pb, i)
#				}
#			}
#			r <- writeStop(r)
#			pbClose(pb)
#			return(r)
#		}		
#	}
#)
#
#
#
#setMethod("Arith", signature(e1='RasterLayerSparse', e2='numeric'),
#    function(e1, e2){ 
#	
#		if (!hasValues(e1)) { stop('RasterLayerSparse has no values') }
#		stopifnot(length(e2) == 1)
#		setValues(e1,  callGeneric(as.numeric(e1@data@values), e2))
#	}
#)
#
#setMethod("Arith", signature(e1='numeric', e2='RasterLayerSparse'),
#    function(e1, e2){ 
#		if (!hasValues(e2)) { stop('RasterLayerSparse has no values') }
#		stopifnot(length(e1) == 1)
#		setValues(e2,  callGeneric(as.numeric(e2@data@values), e1) )
#	}
#)
#
#
#setMethod("Arith", signature(e1='RasterLayer', e2='logical'),
#    function(e1, e2){ 
#		e2 <- as.integer(e2)
#		callGeneric(e1, e2)
#	}
#)
#
#setMethod("Arith", signature(e1='logical', e2='RasterLayer'),
#    function(e1, e2){ 
#		e1 <- as.integer(e1)
#		callGeneric(e1, e2)
#	}
#)
#
#
#
#setMethod("Arith", signature(e1='RasterStackBrick', e2='numeric'),
#    function(e1, e2) {
#	
#		if (length(e2) > 1) {
#			nl <- nlayers(e1)
#			if (length(e2) != nl) {
#				a <- rep(NA, nl)
#				a[] <- e2
#				e2 <- a
#			}
#					
#			if (canProcessInMemory(e1, 3)) {
#				b <- brick(e1, values=FALSE)
#				return( setValues(b, t(callGeneric( t(getValues(e1)), e2))) )
#			}
#			
#			b <- brick(e1, values=FALSE)
#			tr <- blockSize(b)
#			pb <- pbCreate(tr$n)
#			b <- writeStart(b, filename=rasterTmpFile(), bandorder='BIL')
#			for (i in 1:tr$n) {
#				v <- t (callGeneric( t(getValues(e1, row=tr$row[i], nrows=tr$nrows[i])), e2) )
#				b <- writeValues(b, v, tr$row[i])
#				pbStep(pb, i)
#			}
#			b <- writeStop(b)
#			pbClose(pb)
#			return(b)
#		}
#		
#		# else:
#		
#		if (canProcessInMemory(e1, 3)) {
#			b <- brick(e1, values=FALSE)
#			return ( setValues(b,  callGeneric(getValues(e1), e2) ) )
#		} else {
#			b <- brick(e1, values=FALSE)
#			tr <- blockSize(b)
#			pb <- pbCreate(tr$n)
#			b <- writeStart(b, filename=rasterTmpFile())
#			for (i in 1:tr$n) {
#				v <- callGeneric( getValues(e1, row=tr$row[i], nrows=tr$nrows[i]), e2)
#				b <- writeValues(b, v, tr$row[i])
#				pbStep(pb, i)
#			}
#			b <- writeStop(b)
#			pbClose(pb)
#			return(b)
#		}
#	}
#)
#
#
#
#setMethod("Arith", signature(e1='numeric', e2='RasterStackBrick'),
#    function(e1, e2) {
#	
#		if (length(e1) > 1) {
#			nl <- nlayers(e2)
#			if (length(e1) != nl) {
#				a <- rep(NA, nl)
#				a[] <- e1
#				e1 <- a
#			}
#					
#			if (canProcessInMemory(e2, 3)) {
#				b <- brick(e2, values=FALSE)
#				return( setValues(b, t(callGeneric( e1, t(getValues(e2))))) )
#			}
#			
#			b <- brick(e2, values=FALSE)
#			tr <- blockSize(b)
#			pb <- pbCreate(tr$n)
#			b <- writeStart(b, filename=rasterTmpFile())
#			for (i in 1:tr$n) {
#				v <- t (callGeneric( e1, t(getValues(e2, row=tr$row[i], nrows=tr$nrows[i]))) )
#				b <- writeValues(b, v, tr$row[i])
#				pbStep(pb, i)
#			}
#			b <- writeStop(b)
#			pbClose(pb)
#			return(b)
#		}
#		
#		# else:
#		
#		if (canProcessInMemory(e2, 3)) {
#			b <- brick(e2, values=FALSE)
#			return ( setValues(b,  callGeneric(e1, getValues(e2)) ) )
#		} else {
#			
#			b <- brick(e2, values=FALSE)
#			tr <- blockSize(b)
#			pb <- pbCreate(tr$n)
#			b <- writeStart(b, filename=rasterTmpFile())
#			for (i in 1:tr$n) {
#				v <- callGeneric( e1, getValues(e2, row=tr$row[i], nrows=tr$nrows[i]))
#				b <- writeValues(b, v, tr$row[i])
#				pbStep(pb, i)
#			}
#			b <- writeStop(b)
#			pbClose(pb)
#			return(b)
#		}
#	}
#)
#
#
#
#
#
#setMethod("Arith", signature(e1='RasterStackBrick', e2='logical'),  # for Arith with NA
#    function(e1, e2){ 
#		e2 <- as.integer(e2)
#		callGeneric(e1, e2)
#	}
#)
#
#setMethod("Arith", signature(e1='logical', e2='RasterStackBrick'),
#    function(e1, e2){ 
#		e1 <- as.integer(e1)
#		callGeneric(e1, e2)
#	}
#)
#
#
#setMethod("Arith", signature(e1='Extent', e2='numeric'),
#	function(e1, e2){ 
#	
#		if (length(e2) == 1) {
#			x1 = e2
#			x2 = e2
#		} else if (length(e2) == 2) {
#			x1 = e2[1]
#			x2 = e2[2]
#		} else {
#			stop('On an Extent object, you can only use Arith with a single number or with two numbers')
#		}
#
#		r <- e1@xmax - e1@xmin
#		d <- callGeneric(r, x1)
#		d <- (d - r) / 2
#		e1@xmax <- e1@xmax + d
#		e1@xmin <- e1@xmin - d
#		
#		r <- e1@ymax - e1@ymin
#		d <- callGeneric(r, x2)
#		d <- (d - r) / 2
#		e1@ymax <- e1@ymax + d
#		e1@ymin <- e1@ymin - d
#		
#		return(e1)
#	}
#)
#
#setMethod("Arith", signature(e1='numeric', e2='Extent'),
#    function(e1, e2){ 
#		callGeneric(e2,e1)
#	}
#)
#
