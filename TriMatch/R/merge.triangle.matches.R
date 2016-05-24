#' Merges outcomes with the matched set.
#' 
#' The \code{y} parameter should be a subset of the orignal data used.
#' 
#' @param x the result of \code{\link{trimatch}}
#' @param y another data frame or vector to merge with.
#' @param ... unused
#' @return \code{x} with the additional column(s) added.
#' @export
#' @method merge triangle.matches
merge.triangle.matches <- function(x, y, ...) {
	stopifnot('triangle.matches' %in% class(x))
	tpsa <- attr(x, 'triangle.psa')
	
	if(is.vector(y)) {
		if(length(y) != nrow(tpsa)) {
			stop('The length of the vector does not match that of the original data')
		}
		x <- cbind(x, y[as.integer(x[,1])])
		x <- cbind(x, y[as.integer(x[,2])])
		x <- cbind(x, y[as.integer(x[,3])])
		names(x)[(ncol(x)-2):ncol(x)] <- paste(names(x)[1:3], '.out', sep='')
	} else {
		if(nrow(y) != nrow(tpsa)) {
			stop('The length of the vector does not match that of the original data')
		}
		x <- cbind(x, y[as.integer(x[,1]),])
		x <- cbind(x, y[as.integer(x[,2]),])
		x <- cbind(x, y[as.integer(x[,3]),])
		names(x)[(ncol(x)-(3*ncol(y))+1):ncol(x)] <- paste(rep(names(x)[1:3], each=ncol(y)), 
						names(y), sep='.')
	}
	return(x)
}
