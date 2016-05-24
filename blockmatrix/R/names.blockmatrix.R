NULL


#' 
#'\code{names} S3 method for \code{blockmatrix} object
#'
#' @param x a \code{blockmatrix} object 
#' @export
#' @rdname names
#' @method names blockmatrix
#' @S3method names blockmatrix
#' @aliases names
#' @author Emanuele Cordano 
#' 
#' 


names.blockmatrix <- function (x)  {
	
	temp <- x
	class(temp) <- "list"
	out <- names(temp)
	out <- out[out!="value"]
	return(out)
	
}