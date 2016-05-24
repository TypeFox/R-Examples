# TODO: Add comment
# 
# Author: ecor


###############################################################################
NULL

#' @rdname geotopbrick
#' @title geotopbrick
#' @export


geotopbrick <- function(x=NULL,...) {
	
	return(UseMethod("geotopbrick",x))
	
}


NULL

#' 
#' @rdname geotopbrick
#' 
#' @method geotopbrick default
## @S3method geotopbrick default
#' 
#' 
#' @export

geotopbrick.default <- function (x,...) {
	
	out <- new("GeotopRasterBrick")
###	print("default")
	return(out)
	
}


NULL

#'
#' geotopbrick method bla bla bla
#'  
#' @param x a 'zoo' object returned by function \code{\link{pointer.to.maps.xyz.time}} or \code{\link{pointer.to.maps.xy.time}} or a \code{\link{GeotopRasterBrick-class}} object
#' @param layer layer at which raster maps are imported. If is \code{NULL}, maps ara no-zlayer distributed and \code{zoo} must be returend by \code{\link{pointer.to.maps.xy.time}}
#' @param time vector of time instants at which geotop maps are imported
#' @param crs coordinate system see \code{\link{RasterBrick-class}}
#' @param timerange two-elememts vector containing the time range at which geotop maps are imported
#' @param ascpath \code{NULL} object or a \code{"zoo"} S3 object containing the names of ascii maps provided by GEOtop 
#' @param ... further arguments. 
#' 
#' 
#' 
#' @return a \code{\link{GeotopRasterBrick-class}}
#' @title geotopbrick
#' @name geotopbrick
#' 
#' 
#' @rdname geotopbrick
#' 
#' @method geotopbrick zoo
#@S3method geotopbrick zoo
#' @aliases geotopbrick.zoo 
#' 
#' 
#' @export

geotopbrick.zoo <- function(x,layer=NULL,time=NULL,crs=NULL,timerange=NULL,...) { 
		
			if (!is.null(layer)) layer <- 1
			if (is.null(time)) time <- index(x)
			if (!is.null(timerange)) time <- index(x)[index(x)>=timerange[1] & index(x)<=timerange[2]] 
			
			out <- new("GeotopRasterBrick")
			
			out@ascpath <- x
			out@layer <- as.character(layer)
		
			out@brick <- brick(x,layer=layer,time=time,crs=crs,...) 
			out@index <- time		
	
			return(out)
		}


NULL 


#' 
#' @rdname geotopbrick
#' 


#' @method geotopbrick RasterLayer
# @S3method geotopbrick RasterLayer
#' @aliases geotopbrick.RasterLayer
#' 
#' @export

geotopbrick.RasterLayer <- function(x,layer=NULL,time=NULL,ascpath=zoo(NULL),...) {
	
	if (is.null(layer)) layer <- paste("L",1:nlayers(x),sep="")
	if (is.null(time))  time <- as.POSIXlt(character(0))
	
	out <- new("GeotopRasterBrick")
	out@brick <- brick(x,...)
	out@ascpath <- ascpath
	
	
	out@layer <- layer
	
	
	out@index <- time
	
	return(out)
	
}


NULL



#' 
#' @rdname geotopbrick
#' 


#' @method geotopbrick RasterBrick
# @S3method geotopbrick RasterBrick
#' @aliases geotopbrick.RasterBrick
#' 
#' @export

geotopbrick.RasterBrick <- function(x,layer=NULL,time=NULL,ascpath=zoo(NULL),...) {
	
	if (is.null(layer)) layer <- names(x)
	if (is.null(layer)) layer <- paste("L",1:nlayers(x),sep="")
	if (is.null(time))  time <- as.POSIXlt(character(0))
	
	out <- new("GeotopRasterBrick")
	out@brick <- x   ###brick(x,...)
	out@ascpath <- ascpath
	
	
	out@layer <- layer
	
	
	out@index <- time
	
	return(out)
	
}


NULL



#' 
#' @rdname geotopbrick
#' 
#' @method geotopbrick GeotopRasterBrick
# @S3method geotopbrick GeotopRasterBrick
#' 
#' @aliases geotopbrick.GeotopRasterBrick
#' @export
#' 
#' 


geotopbrick.GeotopRasterBrick <- function(x,layer=NULL,time=NULL,crs=NULL,timerange=NULL,ascpath=NULL,...) {
	
	out <- new("GeotopRasterBrick")
	out@brick <- x@brick
	if (is.null(ascpath)) {
		out@ascpath <- x@ascpath
		
	} else {
		
		out@ascpath <- ascpath
		
	}	
	if (is.null(layer)) { 
	
		out@layer <- x@layer
	} else { 
	
		out@layer <- layer
	}
	
	if (is.null(time)) { 
		
		out@time <- x@time
	
	} else {
		
	    out@time <- time
	}
	
	
	return(out)
	
}






