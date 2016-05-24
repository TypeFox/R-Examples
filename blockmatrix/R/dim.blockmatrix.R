NULL

#' 
#'\code{dim} S3 method for \code{blockmatrix} object
#'
#' @param x a \code{blockmatrix} object 
#' @export
#' @rdname dim
#' @method dim blockmatrix
#' @S3method dim blockmatrix
#' @aliases dim
#'
#' @author Emanuele Cordano 
#' 
#' 


dim.blockmatrix <- function (x)  {
	
	return(dim(as.matrix(x$value)))
	
}