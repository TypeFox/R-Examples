# TODO: Add comment
# 
# Author: ecor
###############################################################################


NULL


#'
#' Writes one or more variables (scalars, vectors or Rasters) in a string each, following \code{*.inpts} or Matlab-like syntax. 
#' 
#z-layer brick referred to a time instant (e.g. date) in an ascii format like \code{'geotop.inpts'} file. 
#' 
#' @param l a code{list} object contained the variables (scalars, vectors or Rasters) which will be written in a string each. 
#' @param NAflag numeric. Default is -9999, see \code{\link{writeRasterxGEOtop}}.
#' @param matlab.syntax logical value. Default is \code{FALSE}. If \code{TRUE} the file syntax is like the one of a *.m Matlab script file.
#' @param ... further aguments 
#' @export 
#' 
#' 
#' 
#' @note Add Quote if necessary 
#' @seealso \code{\link{read.ascii.vectorized.brick}}
#' @return the string vector  \code{<NAME_VARIABLE>==<VALUES_VARIABLE>}. 
#' 
#' @examples
#' a <- 1:5
#' l <- list(v=a,a=a)
#' out <- write.vectorized.variable.in.string(l,matlab.syntax=TRUE)
#' out 
#' 
#' 







write.vectorized.variable.in.string <- function(l,NAflag=-9999,matlab.syntax=FALSE,...) {
	

###	l <- list(...)

	if (!is.list(l)) l <- list(l)
	
	nl <- length(l)
	out <- NULL

	if (nl==1) {
	
		x <- l[[1]]
		nv <- names(l)[1]
	
		
		if (length(x)==1) {
			
			x[is.na(x)] <- NAflag
			
			out <- paste(nv,x,sep="=")
			
			return(out)
		} else if (length(x)>1) {
		
			if (class(x)=="RasterBrick" | class(x)=="RasterLayer") {
				vals <- as.vector(getValues(x))
			
			
			} else {
				vals <- x
			}
			
			x <- vals
			
			x[is.na(x)] <- NAflag
			
			if (is.character(x)) x <- paste("\"",x,"\"",sep="")
		
			x <- paste(as.character(x),collapse=",")
			if (matlab.syntax) x <- paste("[",x,"]",sep="")
			out <- paste(nv,x,sep="=")
			
			return(out)
		} 	
		
		return(out)
		
	} else if (nl>1) {
		
		names <- names(l)
		# WARNING: The names of list elements are not saved, so they are 
		out <- lapply(X=l,FUN=write.vectorized.variable.in.string,NAflag=NAflag,matlab.syntax=matlab.syntax,...)
		
		out <- paste(names,unlist(out),sep="")
		
		
		return(out)
	
	}
	
	return(out)
	
}