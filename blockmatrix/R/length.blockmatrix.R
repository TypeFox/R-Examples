NULL

#' 
#'\code{length} S3 method for \code{blockmatrix} object
#'
#' @param x a \code{blockmatrix} object 
#' @export
#' @rdname length
#' @method length blockmatrix
#' @S3method length blockmatrix
#' @aliases length
#'
#' @author Emanuele Cordano 
#' 
#' 


length.blockmatrix <- function (x)  {
	
	return(length(x$value))
	
}