NULL

#' 
#'\code{t} 'transpose' S3 method for \code{blockmatrix} object
#'
#' @param x a \code{blockmatrix} object 
#' @export
#' @rdname t
#' @method t blockmatrix
#' @S3method t blockmatrix
#' @aliases t
#'
#' @author Emanuele Cordano 
#' 
#' 



t.blockmatrix <- function (x)  {

	out <- lapply(FUN=t,x)
	class(out) <- "blockmatrix"
	return(out)
	
}