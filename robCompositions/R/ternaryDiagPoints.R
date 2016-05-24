#' Add points or lines to a given ternary diagram.
#' 
#' Low-level plot function to add points or lines to a ternary high-level plot.
#' 
#' 
#' @aliases ternaryDiagPoints ternaryDiagLines
#' @param x Three-dimensional composition given as an object of class
#' \dQuote{matrix} or \dQuote{data.frame}.
#' @param \dots Additional graphical parameters passed through.
#' @return no values are returned.
#' @author Matthias Templ
#' @seealso \code{\link{ternaryDiag}}
#' @references C. Reimann, P. Filzmoser, R.G. Garrett, and R. Dutter:
#' Statistical Data Analysis Explained. Applied Environmental Statistics with
#' R. John Wiley and Sons, Chichester, 2008.
#' @keywords aplot
#' @export
#' @examples
#' 
#' data(coffee)
#' x <- coffee[,2:4]
#' ternaryDiag(x, grid=FALSE)
#' ternaryDiagPoints(x+1, col="red", pch=2)
#' 
ternaryDiagPoints <- function(x, ...){
	s <- rowSums(x)
	if (any(s <= 0)) 
		stop("rowSums of the input data x must be positive.")
	x <- x/s
	top <- sqrt(3)/2
	xp <- x[, 2] + x[, 3]/2
	yp <- x[, 3] * top
	points(xp, yp, ...)
}
