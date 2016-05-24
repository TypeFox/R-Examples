NULL
#'
#' Extracts daily maximum temparature from an object of class \code{\link{climdexInput-class}}. 
#' 
#' @param x an object of class \code{\link{climdexInput-class}}
#' @return Daily Maximum Temperature 
#' 
#' @author Emanuele Cordano, Annalisa Di Piazza
#' @title Daily Maximum Temperature
#' @seealso \code{\link{climdexInput-class}}, \code{\link{climdexInput.raw}}
#' @export

temperature_max_daily <- function(x) {
	
	print(x@tmax)
	return(x@tmax)
	
	
}