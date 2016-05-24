# TODO: Add comment
# 
# Author: ecor
###############################################################################


NULL


#'
#' Read a text file containing  values and matedata of a z-layer brick referred to a time instant (e.g. date). The file is formatted like an ascii format like \code{'geotop.inpts'} file. 
#' 
# @param b a \code{\link{RasterBrick-class}} or \code{\link{GeotopRasterBrick-class}} object
#' @param file file name to write 
#' @param comment character. Comment indicator. Default is \code{"!"}.
# @param header character string vector for header text lines. Default is  \code{c("! header")}. 
# @param overwrite logical. Default is \code{TRUE}, see \code{\link{writeRaster}}.
#' @param NAflag numeric. Dafauli is -9999, see \code{\link{writeRasterxGEOtop}}.
#' @param crs Character or object of class CRS. PROJ4 type description of a Coordinate Reference System (map projection) (optional). See \code{\link{brick}} or \code{\link{raster}}.
#' @param matlab.syntax logical value. Default is \code{FALSE}. If \code{TRUE} the file syntax is like the one of a *.m Matlab script file.
#' @param ... further aguments inserted as attribute
# @export 
#' 
#' 
#' 
#' @return the \code{\link{RasterBrick-class}} object 
#' @export
#' @seealso \code{\link{write.ascii.vectorized.brick}}
#' 
#' @examples 
#' # see the examples of read.ascii.vectorized.brick





read.ascii.vectorized.brick <- function(file=NULL,comment="!",crs="",NAflag=-9999,matlab.syntax=FALSE,...) {
	
	if (matlab.syntax) comment="#"

	df <- declared.geotop.inpts.keywords(inpts.file=file,comment=comment,warn=FALSE,wpath=NULL,...) # commented wpath='.'
	
	xmx <- get.geotop.inpts.keyword.value("xmx",inpts.frame=df,numeric=TRUE)
	xmn <- get.geotop.inpts.keyword.value("xmn",inpts.frame=df,numeric=TRUE)
	ymx <- get.geotop.inpts.keyword.value("ymx",inpts.frame=df,numeric=TRUE)
	ymn <- get.geotop.inpts.keyword.value("ymn",inpts.frame=df,numeric=TRUE)
	
	nrow <- get.geotop.inpts.keyword.value("nrow",inpts.frame=df,numeric=TRUE)
	ncol <- get.geotop.inpts.keyword.value("ncol",inpts.frame=df,numeric=TRUE)
	nlayers <- get.geotop.inpts.keyword.value("nlayers",inpts.frame=df,numeric=TRUE)
	## IF MATLAB 
	brickvalues <- get.geotop.inpts.keyword.value("brickvalues",inpts.frame=df,numeric=TRUE)
	
	nh <- nrow*ncol*nlayers
	
	if (length(brickvalues)<nh) {
		
		print("Error in read.ascii.vectorized.formatter: ncol or nrow or nlayers are incompatible with brickvalues!!")
		stop()
		
	}
	
	b <- brick(nrows=nrow, ncols=ncol, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=crs, nl=nlayers)
# convention in x-y-z 
#	x <<- array(NA,c(nrow,ncol,nlayers))
#	print(class(x))
	for (i in 1:nlayers) {
		
		element <- nrow*ncol*(i-1)+1:(nrow*ncol)
		
		M <- (array(brickvalues[element],c(ncol,nrow)))
		b <- setValues(b,M,layer=i)
	#	print(M)
	#	print(x[1:nrow,1:ncol,i])
	#	x[1:nrow,1:ncol,i] <- M[1:nrow,1:ncol]
	#	b[[i]][,] <- M
	#	value(b[[i]]) <- M
		
	}
	
#	b <- brick(x=brickvalues, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=crs)
#	print(xmx)
#	print(xmn)
#	print(ymx)
#	print(ymn)
	
	
	
	out <- b
#	l <- list(...)
	
#	if (class(b)=="GeotopRasterBrick") b <- brick(b)
#	
#	l$xmx <- xmax(b)
#	l$ymx <- ymax(b)
#	l$xmn <- xmin(b)
#	l$ymn <- ymin(b)
#	l$nlayers <- nlayers(b)
#	l$nrow <- nrow(b)
#	l$ncol <- ncol(b)
##  X Y Z
#	vals <- as.vector(getValues(b))
#	vals[is.na(vals)] <- NAflag
#	l$brickvalues <- paste(as.character(vals),collapse=",")
#	
##	lines <- c(header,array("",length(l)))
##	lines <- array("",length(l))
#	
#	lines <- names(l)
#	arg <- unlist(l)
#	
#	lines <- paste(lines,arg,sep="=")
##	for ( i in 1:length(lines)) {
#		
##		lines[i] <- paste
#		
##	}
#	
#	
#	
#	lines <- c(header,lines)
#	
#	if (!is.null(file)) {
#		
#		writeLines(lines, con = file, sep = "\n", useBytes = FALSE)
#	}
	
	return(out)
	
}