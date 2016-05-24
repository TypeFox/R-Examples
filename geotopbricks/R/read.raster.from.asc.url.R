NULL
#'
#' It imports a 'RasterLayer' object in Escri-Asci format from a URL 'http://....<FILENAME>.asc
#' 
#' @param x the charcater string containing the URL address
#' @param header_nrow Number of header in the ASCII grid format. Deafault is 6. See \url{http://en.wikipedia.org/wiki/Esri_grid}
#' @param ... additional arguments
#'
#' @note This function reads a local or remote text files formatted as \url{http://en.wikipedia.org/wiki/Esri_grid} and creates a 'RasterLayer' object.
#' @export 
#' @return a 'RasterLayer' object
#' @seealso \code{\link{raster}},\code{\link{readLines}}
###,\url{http://en.wikipedia.org/wiki/Esri_grid}
#' 

read.raster.from.url <- function(x,header_nrow=6,...) {
	
#	out <- NULL
	
#	y <- url(x,...)
	
#	# ec date 09-03-2013
#	if (str_sub(x,1,3)=='ssh' | str_sub(x,1,5)=='plink') {
#		
#		file <- pipe(x) # added line according to http://stackoverflow.com/posts/2226880/edit
#		 open <- TRUE                   
#	}	else {
#		
#		file <- file(x)
#		open <- FALSE
#	}
#	# ec date 09-03-2013
###	file <- file(x)    ### 2014-05-20 ec 
	lin <- readLines(x,warn=FALSE)
	
##	if (open) close(file) # added ec date 09-03-2013
	
	s0 <- str_split(lin,pattern=" ")
	
	header <- unlist(s0[1:header_nrow])
	header <- header[header!=""]
	
	ncols <- as.numeric(header[which(header=="ncols" | header=="NCOLS")+1])
	nrows <- as.numeric(header[which(header=="nrows" | header=="NROWS")+1])
	cellsize <- as.numeric(header[which(header=="cellsize" | header=="CELLSIZE")+1])
	xllcorner <- as.numeric(header[which(header=="xllcorner" | header=="XLLCORNER")+1])
	yllcorner <- as.numeric(header[which(header=="yllcorner" | header=="YLLCORNER")+1])
	NA_data <- as.numeric(header[which(header=="NODATA_value")+1])
	

	values0 <- lapply(FUN=as.numeric,s0[-c(1:header_nrow)])
	values <- array(unlist(values0),c(ncols,nrows))
	values <- t(values)

	values[values==NA_data] <- NA 

	xmn <- xllcorner
	ymn <- yllcorner
	
	xmx <- xmn+cellsize*ncols
	ymx <- ymn+cellsize*nrows
	
	out <- raster(x=values,xmn=xmn,ymn=ymn,xmx=xmx,ymx=ymx,...)

	
	return(out)
	
}