#' Merges covariate(s) with the results of \code{\link{trips}}.
#' 
#' The \code{y} parameter should be a subset of the orignal data used.
#' 
#' @param x the result of \code{\link{trips}}
#' @param y another data frame or vector to merge with.
#' @param ... unused
#' @return \code{x} with the additional column(s) added.
#' @export
#' @method merge triangle.psa
merge.triangle.psa <- function(x, y, ...) {
	stopifnot('triangle.psa' %in% class(x))
	if(is.vector(y) & length(y) != nrow(x)) {
		stop('The length of the vector does not match that of the original data')
	} else if(nrow(y) != nrow(x)) {
		stop('The length of the vector does not match that of the original data')
	}
	x <- cbind(x, y)
	return(x)
}
