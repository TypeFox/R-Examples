#' Method for creating xyloplots
#'
#' @param x Numeric vector or list of numeric vectors to use for creating xyloplots.
#' @param ... Other arguments to be passed to \code{\link{xyloplot.list}} and \code{plot}.
#' @examples
#' xyloplot(rnorm(1000))
#' xyloplot(
#'  x=lapply(1:3, function(mean) rnorm(mean=mean, n=1000)), 
#'  breaks=20,
#'  col=rainbow(3), 
#'  main="title")
#' @seealso xyloplot.list xyloplot.numeric
#' @export
xyloplot <- function(x, ...) {
	UseMethod("xyloplot")	
}

#' Function for creating xyloplots from numeric vectors
#'
#' @param x Numeric vector of values.
#' @param ... Other arguments passed to \code{\link{xyloplot.list}}.
#'
#' @method xyloplot numeric
#' @seealso xyloplot.list
#' @export
xyloplot.numeric <- function(x, ...) {
	xyloplot.list(list(x), ...)
}

#' Function for creating multiple xyloplots sharing a y-axis from lists of numeric vectors
#'
#' @param x List of numeric vectors of values.
#' @param ylim Limits of hte y-axis.
#' @param breaks A single positive integer value giving the number of histogram classes to evenly split the values in \code{x} into, or a numeric vector explicitly giving the boundaries of the histogram classes.
#' @param space The proportion of the total distance on the x-axis allocated to each 'xylophone' which should be left blank.
#' @param ylab Label for y-axis.
#' @param xlab Label for x-axis.
#' @param ... Other arguments to be passed to \code{\link{plot}}.
#' @method xyloplot list
#' @export
#' @importFrom graphics axis box plot rect title
xyloplot.list <- function(x, ylim=range(unlist(x, use.names=FALSE)), breaks=30, space=0.1, ylab="Value", xlab="Frequency density", ...) { 
	stopifnot(length(breaks) > 1 | (length(breaks) == 1 & breaks > 1))
	stopifnot(is.numeric(unlist(use.names=FALSE, x)))
	stopifnot(diff(range(unlist(use.names=FALSE, x)))>0)

	brk.pts <- if (length(breaks) == 1) seq(from=ylim[1], to=ylim[2], length.out=breaks) else breaks
	
	blks <- lapply(FUN=function(x) as.array(table(x)/length(x)), lapply(x, cut, breaks=brk.pts, right=FALSE))

	plot(x=NULL, xlim=0:1, xaxs="i", yaxs="i", xlab="", ylab=ylab, ylim=if (is.null(ylim)) range(brk.pts) else ylim, axes=FALSE, ...)

	violin.width <- (1-space)/length(x)

	pivots <- seq(from=1/2/length(x),by=1/length(x),length.out=length(x))
	ord <- as.integer(t(matrix(1:(length(x) * (length(brk.pts)-1)), nrow=length(brk.pts)-1,ncol=length(x))))

	rect(
		xleft=(rep(pivots, each=length(brk.pts)-1)-unlist(blks)/max(unlist(blks))*violin.width/2)[ord],
		xright=(rep(pivots, each=length(brk.pts)-1)+unlist(blks)/max(unlist(blks))*violin.width/2)[ord],
		ybottom=(rep(brk.pts[-1], times=length(x))-rep(diff(brk.pts), times=length(x)))[ord],
		ytop=(rep(brk.pts[-1], times=length(x)))[ord],
		...
	)

	box()

	if (!is.null(names(x))) {
		axis(tick=FALSE, side=1,at=pivots,labels=names(x))
		title(xlab=xlab, line=3)
	} else {
		title(xlab=xlab, line=1)
	}

	axis(side=2)
}
