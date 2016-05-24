# TODO: Add comment
# 
# Author: ecor
###############################################################################


NULL


#'
#' It writes a \code{list} object returened by \code{\link{get.geotop.recovery.state}} as a string vector or in a text file, following \code{*.inpts} or Matlab-like syntax.  
#' 
#' @param rec a \code{list} object returened by \code{\link{get.geotop.recovery.state}}
#' @param file ascii text file name whrere to write the string vector
#' @param header character string vector for header text lines. If missing, a default header is written. Default is  \code{c("! header")} or he one assigned by \code{matlab.syntax}.
#' @param overwrite logical. Default is \code{TRUE}, see \code{\link{writeRaster}}.
#' @param NAflag numeric. Default is -9999, see \code{\link{writeRasterxGEOtop}}.
#' @param matlab.syntax logical value. Default is \code{TRUE}. If \code{TRUE} the file syntax is like the one of a *.m Matlab script file.
#' @param ... further aguments inserted as attribute
#' @export
#' 
#' @note Add Quote if necessary 
#' @seealso \code{\link{get.geotop.recovery.state}},\code{\link{set.geotop.recovery.state}},\code{\link{write.vectorized.variable.in.string}}
#' @return a string vector containg the \code{rec} variables. 
#' 
#' @examples 
#' # See the examples of the 'get.geotop.recovery.state' function
#' 
#' 
#' 









write.vectorized.geotop.recovery <- function(rec,file=NULL,header=NULL,overwrite=TRUE,NAflag=-9999,matlab.syntax=TRUE,...) {
	

	is.rec <- (class(rec)=="list")
	
	rec_attributes <- c("names","files","noLayers","soilLayersWithZero","soilLayers","snowLayers") 
	is.rec <- is.rec & as.logical(min(rec_attributes %in% names(rec),na.rm=TRUE))
##	is.rec <- is.rec & !(is.null(rec$names) | is.null(rec$files) | is.null(rec$noLayers) | is.null(rec$soilLayersWithZero) | is.null(rec$soilLayers) | is.null(rec$snowLayers)) 
	
#	if (is.rec) {
#		
#		names <- is.rec$names
#		for (it in names) {
#			if()
#		}
#	}
	
	
	if (!is.rec) {

		print("Error in write.vectorized.geotop.recovery: rec is not a recovery list of GEOtop!!!")
		stop()
	}
	names <- rec$names
	rec_attributes0 <- rec_attributes[!(rec_attributes %in% c("names","files"))]
	l <- list(...)

	if (is.null(header)) {
		
		f <- system.file("doc/examples/snowthickness_textfile/template.txt",package="geotopbricks")
		if (matlab.syntax) f <- system.file("doc/examples/snowthickness_textfile/template_matlab.txt",package="geotopbricks")
		header <- readLines(f,warn=FALSE)
	}
	
	
	postheader <- paste("## DISTRUBUTED RASTER VARIABLES: ")
	postheader <- array(postheader,length(rec_attributes0)+1)
	smh <- 6
	metadataheader <- array(" ",length(rec_attributes0)+smh)
	
	metadataheader[1:smh] <- "THIS IS TO DO!!!!"
	
	
	xmx <- xmax(rec[[names[1]]])
	xmn <- xmin(rec[[names[1]]])
	
	ymx <- ymax(rec[[names[1]]])
	ymn <- ymin(rec[[names[1]]])
	
	nrow <- nrow(rec[[names[1]]])
	ncol <- ncol(rec[[names[1]]])
	
	metadataheader[1] <- paste("xmx=",xmx,sep="")
	metadataheader[2] <- paste("xmn=",xmn,sep="")
	metadataheader[3] <- paste("ymx=",ymx,sep="")
	metadataheader[4] <- paste("ymn=",ymn,sep="")
	
	metadataheader[5] <- paste("nrow=",nrow,sep="")
	metadataheader[6] <- paste("ncol=",ncol,sep="")
	###xmx <- xmax(b)
	for (i in 1:length(rec_attributes0)) {
	
		n <- rec_attributes0[i]
		v <- rec[[n]]
		nms <- names[v]
###		postheader[i+1] <- paste("#!!",n,":",paste(nms,collapse=","),sep=" ")
##		postheader[i+1] <- paste(n,paste(nms,collapse=","),sep="=") #RIMETTERE A POSTO QUI !!!!
		li <- list(x=nms)
		names(li) <- n
		
###		print(li)
		postheader[i+1] <- write.vectorized.variable.in.string(li,NAflag=NAflag,matlab.syntax=matlab.syntax,...)
		# TEST THIS FUNCTION!!!!
		
		b <- rec[[nms[1]]]
		
		metadataheader[i+smh] <- paste("nlayers_",n,"=",nlayers(b),sep="")
		
	}
	
	
	l[names] <- rec[names]
	
	

	
	lines <- write.vectorized.variable.in.string(l,NAflag=NAflag,matlab.syntax=matlab.syntax,...)
	
	# ADD METADATA
	
	lines <- c(header,postheader,metadataheader,lines)
	



#	l <- list(...)
#	
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
#	if (matlab.syntax) l$brickvalues <- paste("[",l$brickvalues,"]",sep="")
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
#	if (is.null(header)) {
#		
#		f <- system.file("doc/examples/snowthickness_textfile/template.txt",package="geotopbricks")
#		if (matlab.syntax) f <- system.file("doc/examples/snowthickness_textfile/template_matlab.txt",package="geotopbricks")
#		header <- readLines(f,warn=FALSE)
#	}
#	
#	lines <- c(header,lines)
#	
	if (!is.null(file)) {
		
		writeLines(lines, con = file, sep = "\n", useBytes = FALSE)
	}
#	
#	return(lines)
	
	
	
	
}