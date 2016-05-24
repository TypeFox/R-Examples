NULL 
#'
#' 
#' 

#' @rdname value
#' @export
value <- function (M)  {
	
	return(UseMethod("value",M))
	
}


NULL
#' 
#' 

#' 
#' @rdname value
#' @method value default
#' @S3method value default
#' @aliases value 
#' @export



value.default <- function (M)  {
	
	return(NULL)
	
}
NULL

#' 
#'\code{value} S3 method for \code{blockmatrix} object
#'
#' @param M a \code{blockmatrix} object 
#' 
#' @return The character matrix without numerical values (e.g. only the matrix \code{M$value})
#' 
#' @export
#' @rdname value
#' @method value blockmatrix
#' @S3method value blockmatrix
#' @aliases value
#' @author Emanuele Cordano 
#' 
#' 


value.blockmatrix <- function (M)  {
	
	
	return((M$value))
	
}