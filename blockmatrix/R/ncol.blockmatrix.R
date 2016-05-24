NULL


#' 
#'\code{ncol} S3 method for \code{blockmatrix} object
#'
#' @param M a \code{blockmatrix} object 
#' 
#' 
#' @return Numbner of columns of blockmatrix \code{M} 
#' 
#' @export
#' @rdname ncol
#' @method ncol blockmatrix
#' @S3method ncol blockmatrix
#' @aliases ncol
#' @author Emanuele Cordano 
#' 
#' 


ncol.blockmatrix <- function (M)  {
	
	return(ncol(as.matrix(M$value)))
	
}