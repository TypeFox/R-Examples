#' Paint window for label ?
#' 
#' Paint window for label based on projected coordinates?.
#' 
#' 
#' @param listi Label location ?
#' @param col Color, not used ?
#' @param border Should border be drawn?
#' @param poly Should label be opaque ?
#' @return No value, lines and/or polygon for label added to current geoplot.
#' @note Needs elaboration and possibly merging with paint.window doc-file.
#' Argument \code{col} has no effect.
#' @seealso Called by \code{\link{geocontour}}.
#' @keywords aplot
#' @export paint.window.x
paint.window.x <-
function(listi, col = 0., border = T, poly = T)
{
	x <- list(y = c(listi$y[1.], listi$y[1.], listi$y[2.], listi$y[2.],
		listi$y[1.]), x = c(listi$x[1.], listi$x[2.], listi$x[2.],
		listi$x[1.], listi$x[1.]))
	rx <- range(x$x)
	ry <- range(x$y)
	t1 <- c(rx[1.], rx[2.], rx[2.], rx[1.], rx[1.])
	t2 <- c(ry[1.], ry[1.], ry[2.], ry[2.], ry[1.])
	if(border) {
		mx <- mean(t1[1.:4.])
		my <- mean(t2[1.:4.])
		t11 <- t1 + 0.02 * (t1 - mx)
		t22 <- t2 + 0.02 * (t2 - my)
		lines(t11, t22, lwd = 1.5, col = 1.)
	}
	if(poly)
		polygon(t1, t2, col = 0.)
}

