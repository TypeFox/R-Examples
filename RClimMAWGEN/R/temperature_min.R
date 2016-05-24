NULL
#'
#' Extracts daily Minimum temparature from an object of class \code{\link{climdexInput-class}}. 
#' 
#' @param x an object of class \code{\link{climdexInput-class}}
#' @return Daily Minimum Temperature 
#' 
#' @author Emanuele Cordano, Annalisa Di Piazza
#'
#' @title Daily Minimum Temperature
#'  
#' @seealso \code{\link{climdexInput-class}}, \code{\link{climdexInput.raw}}
#' @export


temperature_min_daily <- function(x) {
	
	return(x@tmin)
	
	
}