NULL

#' 
#'\code{nrow} S3 method for \code{blockmatrix} object
#'
#' @param M a \code{blockmatrix} object 
#' @export
#' @rdname nrow
#' @method nrow blockmatrix
#' @S3method nrow blockmatrix
#' @aliases nrow
#'
#' @return Number of rows of blockmatrix \code{M} 
#' 
#' @author Emanuele Cordano 
#' 
#' 



nrow.blockmatrix <- function (M)  {
	
	return(nrow(as.matrix(M$value)))
	
}