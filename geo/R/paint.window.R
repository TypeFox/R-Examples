#' Paint window for label?
#' 
#' Paint window for label based on lat and lon?.
#' 
#' 
#' @param listi Label location ?
#' @param col Color, not used ?
#' @param border Should border be drawn?
#' @param poly Should label be opaque ?
#' @param col.names Column names containing label coordinates ?
#' @return No value, lines and/or polygon added to current geoplot.
#' @note Needs elaboration and possibly merging with paint.window doc-file.
#' Argument \code{col} has no effect.
#' @seealso Called by \code{\link{colsymbol}}, \code{\link{geocontour}},
#' \code{\link{geocontour.fill}}, \code{\link{geosymbols}} and
#' \code{\link{reitaplott}}; calls \code{\link{Proj}}.
#' @keywords aplot
#' @export paint.window
paint.window <-
function(listi, col = 0, border = T, poly = T, col.names = c("lon", "lat"))
{
	geopar <- getOption("geopar")
	lat <- c(listi[[col.names[2]]][1], listi[[col.names[2]]][1], listi[[
		col.names[2]]][2], listi[[col.names[2]]][2], listi[[col.names[
		2]]][1])
	lon <- c(listi[[col.names[1]]][1], listi[[col.names[1]]][2], listi[[
		col.names[1]]][2], listi[[col.names[1]]][1], listi[[col.names[
		1]]][1])
	x <- Proj(lat, lon, geopar$scale, geopar$b0, geopar$b1, geopar$l1,
		geopar$projection, col.names = col.names)
	rx <- range(x$x)
	ry <- range(x$y)
	t1 <- c(rx[1], rx[2], rx[2], rx[1], rx[1])
	t2 <- c(ry[1], ry[1], ry[2], ry[2], ry[1])
	if(poly)
		polygon(t1, t2, col = 0)
	if(border) {
		mx <- mean(t1[1:4])
		my <- mean(t2[1:4])
		t11 <- t1 + 0.02 * (t1 - mx)
		t22 <- t2 + 0.02 * (t2 - my)
		lines(t11, t22, lwd = 1.5, col = 1)
	}
}

