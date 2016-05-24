# TODO: Add comment
# 
# Author: ecor
###############################################################################


NULL


#'
#' Writes a z-layer brick referred to a time instant (e.g. date) in an ascii format like \code{'geotop.inpts'} file. 
#' 
#' @param b a \code{\link{RasterBrick-class}} or \code{\link{GeotopRasterBrick-class}} object
#' @param file file name to write 
#' @param header character string vector for header text lines. If missing, a default header is written. #Default is  \code{c("! header")}. 
#' @param overwrite logical. Default is \code{TRUE}, see \code{\link{writeRaster}}.
#' @param NAflag numeric. Default is -9999, see \code{\link{writeRasterxGEOtop}}.
#' @param matlab.syntax logical value. Default is \code{FALSE}. If \code{TRUE} the file syntax is like the one of a *.m Matlab script file.
#' @param ... further aguments inserted as attribute
#' @export 
#' 
#' 
#' 
#' @note Add Quote if necessary. This function is NOT mantained and will be DEPRECATED.
#'  
#' @seealso \code{\link{read.ascii.vectorized.brick}}
#' @return the string vector possibly written in \code{file}. 
#' 
#' @examples 
#' ## Not Run
#' ## library(geotopbricks)
#' ## library(raster)
# b <- brick(system.file("external/rlogo.grd", package="raster"))
#' ## file <- system.file("doc/examples/snowthickness",package="geotopbricks")
#' ## file <- paste(file,"SnowThickness0000L%04d.asc",sep="/")
#' ## b <- brick.decimal.formatter(file=file,nlayers=15)
#' ## nlayers(b)
#' ## names(b)
#' ## file <- "snow.txt"
#' ## btext <- write.ascii.vectorized.brick(b,Date="1/1/2009",file="snow.txt")

#' ## The printed object
#' ## str(btext)
#' ## bb <- read.ascii.vectorized.brick(file = file) 
#' ## bf <- abs(as.matrix(bb[[1]]-b[[1]]))<.Machine$double.eps^0.5








write.ascii.vectorized.brick <- function(b,file=NULL,header=NULL,overwrite=TRUE,NAflag=-9999,matlab.syntax=FALSE,...) {
	
	
	l <- list(...)
	
	if (class(b)=="GeotopRasterBrick") b <- brick(b)

	l$xmx <- xmax(b)
	l$ymx <- ymax(b)
	l$xmn <- xmin(b)
	l$ymn <- ymin(b)
	l$nlayers <- nlayers(b)
	l$nrow <- nrow(b)
	l$ncol <- ncol(b)
#  X Y Z
	vals <- as.vector(getValues(b))
	vals[is.na(vals)] <- NAflag
	l$brickvalues <- paste(as.character(vals),collapse=",")
	if (matlab.syntax) l$brickvalues <- paste("[",l$brickvalues,"]",sep="")
	
#	lines <- c(header,array("",length(l)))
#	lines <- array("",length(l))
	
	lines <- names(l)
	arg <- unlist(l)
	
	lines <- paste(lines,arg,sep="=")
#	for ( i in 1:length(lines)) {
		
#		lines[i] <- paste
		
#	}
	
	if (is.null(header)) {
		
		f <- system.file("doc/examples/snowthickness_textfile/template.txt",package="geotopbricks")
		if (matlab.syntax) f <- system.file("doc/examples/snowthickness_textfile/template_matlab.txt",package="geotopbricks")
		header <- readLines(f,warn=FALSE)
	}
	
	lines <- c(header,lines)
	
	if (!is.null(file)) {
		
		writeLines(lines, con = file, sep = "\n", useBytes = FALSE)
	}
	
	return(lines)
	
}