NULL
#'
#' Gets the minimum (scalar) values of a \code{\link{GeotopRasterBrick}} object 
#' 
#' @param x a \code{\link{GeotopRasterBrick}} object 
#' @param ... further arguments 
#' 
#' 
#' @return the minimum (scalar) values of a \code{\link{GeotopRasterBrick}} object
#' 
#' 
#' 
#' @title min_value
#' @name min_value
#' @rdname min_value
#' @export 
#'

min_value <- function(x) {
	
	out <- min(minValue(brick(x)),na.rm=TRUE)
	
	return (out) 
	
}


