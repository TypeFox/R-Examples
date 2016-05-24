NULL 
#'
#' is.zero.bolockmatrix
#' 
#' @param M a \code{blockmatrix} object
#' @param not.a.blockmatrix value to be returned in case \code{M} is not a a \code{blockmatrix} object
#' 
#' @return logical value in case \code{M} is a zero blockmatrix
#' @author Emanuele Cordano
#' @export
is.zero.blockmatrix <- function (M,not.a.blockmatrix=FALSE)  {
	
	if (class(M)!="blockmatrix") return(not.a.blockmatrix)
	v <- value(M)
	out <- FALSE
	if (length(v)==1) {
		if (v==0) out <- TRUE
	}
	
	return(out)
	
}
