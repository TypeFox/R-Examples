NULL

#' Unit Offsets of a \code{\link{smet}} object. 
#' 
#'@param x a \code{\link{smet-class}} object 
#' 
#' @export 
#' 
#' @examples
#' 
#' x <- smet(system.file('examples/PIEM001114.smet',package="RSMET"))
#' units_offset(x)
#' 

units_offset <- function(x) {
	
	return(x@header$units_offset)
	
}






