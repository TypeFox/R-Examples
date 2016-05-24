#' @export
#' @rdname plots
plotcv <- function(obj, stack = FALSE, index.res = 1, col = index.res, ...) { 
	if (stack)
		stop("Unimplemented option: stack!");
	if (index.res != 1)
		stop ("Unimplemented option: index.res");
	if (class(obj)!="LHS")
		stop("The first argument should be of class LHS!");
	if (get.repetitions(obj)<2)
		stop("Error in function plotcv: the LHS object must have at least two repetitions!")
	pointwise <- apply(get.results(obj, FALSE), c(1,2), cv)
	global <- cv(get.results(obj, TRUE))
	m <- max(pointwise, 1.05*global)
	mi <- min(pointwise, global)
	Ecdf(pointwise, xlim=c(mi, m), xlab="pointwise cv", col=col, ...)
	abline(v=global, lwd=2, lty=3)
	if (m > 0.8*max(pointwise)) {pos=2} else {pos=4}
	text(x=global, y=0.1, label="global cv", pos=pos)
}

#' Coefficient of Variation
#' 
#' Returns the coefficient of variation of a sample.
#' @param x Any numeric vector (or other data type for which \code{sd} and \code{mean} methods are defined)
#' @param \dots Additional parameters for the \code{sd} and \code{mean} functions (such as \code{na.rm=TRUE})
#' @export
cv <- function(x, ...) {
	if (mean(x, ...) == 0) return (0)
	else return (sd(x, ...)/mean(x, ...))
}
