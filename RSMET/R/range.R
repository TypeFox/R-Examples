NULL 
#' 'range' method for a \code{\link{smet-class}} object
#' 
#' @param x a \code{\link{smet-class}} object
#' @param field field (character string) used for range. Default is \code{"timestamp"}.
#' @param ... further arguments
#' 
#' 
#' @export
#' 
#' @rdname range
#' @method range smet
#' @aliases range  
#' 

range.smet <- function(x,...,field="timestamp") {
	
	out <- base::range(x@data[,field],...)
	
	return(out)
	
	
	
}