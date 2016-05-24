NULL
#'
#' @rdname value.replace
#' 
#' 
#' @export
#' @usage value(M) <- value



'value<-' <- function (M,value)  {
	
	
	return(UseMethod("value<-",M))
}

NULL
#'
#' @rdname value.replace
#' @method value<- default
#' @S3method value<- default
#' @aliases value<-
#' 
#' @usage \method{value}{default} (M) <- value 
#' @export



'value<-.default' <- function (M,value)  {
	
	
	return(value)
}





NULL
#' 
#'\code{value<-} S3 Replacement method for \code{blockmatrix} object
#'
#' @param M a \code{blockmatrix} object 
#' @param value object replaced matrix 
#' @export
#' 
#' @return Replaces \code{M$value} with a new matrix \code{value}
#' 
#' @rdname value.replace
#' @method value<- blockmatrix
#' @S3method value<- blockmatrix
#' @aliases value<-
#' @author Emanuele Cordano 
#' @usage \method{value}{blockmatrix} (M) <- value
#' 



'value<-.blockmatrix' <- function (M,value)  {
	
	M$value <- value
	return(M)
}