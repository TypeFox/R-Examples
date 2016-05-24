NULL
#'
#' Gets the maximum (scalar) values of a \code{\link{GeotopRasterBrick}} object 
#' 
#' @param x a \code{\link{GeotopRasterBrick}} object 
#' @param ... further arguments 
#' 
#' 
#' @return the maximum (scalar) values of a \code{\link{GeotopRasterBrick}} object
#' 
#' 
#' 
#' @title max_value
#' @name max_value
#' @rdname max_value
#' @export 
#'

max_value <- function(x) {
	
	out <- max(maxValue(brick(x)),na.rm=TRUE)
	
	return (out) 
	
}


