#' Clear Memory of All Objects
#'
#' This function is a wrapper for the command \code{rm(list=ls())}.
#'
#' @param obj The object (as a string) that needs to be removed (or kept)
#' @param keep Should \code{obj} be kept (i.e., everything but \code{obj} removed)? Or dropped?
#' @author Daniel Marcelino
#' @export
#' @examples
#' # create objects
#' a=1; b=2; c=3; d=4; e=5
#' # remove d
#' clear("d", keep=FALSE)
#' ls()
#' # remove all but a and b
#' clear(c("a", "b"), keep=TRUE)
#' ls()
`clear` = function(obj=NULL, keep=TRUE){
	if (!is.null(obj)){
		if (keep){
			dropme = ls(envir=globalenv())[which(!(ls(envir=globalenv())%in%obj))]
		} else {
			dropme = obj
		}
		rm(list=dropme, envir=globalenv())
		cat("All objects were deleted, except:", dropme, sep=",")
	} else {
		rm(list=ls(envir=globalenv()), envir=globalenv())
		cat("All objects were deleted, including hidden package environments.\n")
	}
}


