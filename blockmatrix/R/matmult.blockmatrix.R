NULL

#' 
#'\code{blockmatmult} implements the implents betwwen two blockmatrix ( see \code{\link{matmult}} for \code{matrx} objects)
#' 
#' @name blockmatmult
#' @param x,y  \code{blockmatrix} objects 
#' @param ... further arguments 
#' @export
#' @rdname blockmatmult
# @method matmult blockmatrix
# @S3method matmult blockmatrix
#' @aliases blockmatmult
#'
#' @return The inner product between \code{x} and \code{y} as a \code{blockmatrix} object
#' 
#' @author Emanuele Cordano 
#' 
#' 



blockmatmult  <- function (x,y,...)  {
	
	ncol <- ncol(y)
	nrow <- nrow(x)
	
	xm <- as.matrix(x,...)
	ym <- as.matrix(y,...)
	
	
	
	out <- xm %*% ym 
	
	
	return(as.blockmatrix(out,ncol=ncol,nrow=nrow,...))
	
}