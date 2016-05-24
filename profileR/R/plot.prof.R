#' Plots a pattern and level reliability
#' 
#' Plots the pattern vs. level reliability returned from the \code{pr} function of class \code{prof}.
#' @importFrom graphics plot
#' @param x an object returned from the \code{pr} function
#' @param ... additional objects of the same type.
#' @importFrom graphics plot
#' @method plot prof
#' @seealso \code{\link{pr}}
#' @export
#' 
plot.prof <- function(x, ...){
	dat0 <- x$pattern.level
	dim <- ncol(dat0)
	dat1 <- as.data.frame(dat0[,c(dim/2,dim)])
	plot(dat1, xlab="Level 1", ylab="Level 2", ...)
}


