#' Fills space outside the border of a plot.
#' 
#' When programs like geopolygon are used they sometimes fills space outside
#' the border of their plots, to refill that space white we use
#' fill.outside.border.
#' 
#' 
#' @param col The color of the fill, default is 0 (usually white).
#' @param rat Ratio between the size of the current plot and the size of the
#' fill, if you want to allow the program to draw 5\% outside the plot you set
#' rat =1.05.  Default rat = 1.
#' @return No Value.
#' @section Side Effects: Fill outside of the border of the current plot.
#' @seealso \code{\link{geoworld}}, \code{\link{geopolygon}}.
#' @keywords <!--Put one or more s-keyword tags here-->
#' @examples
#' 
#'        \dontrun{
#' 	geoplot(xlim=c(-50,20),ylim=c(50,70))       # Initialize plot.
#'        geoworld(fill=T,color=120)                  # Colour countries.
#'        fill.outside.border()                       # Clear outside of border.
#'        geoplot(xlim=c(-50,20),ylim=c(50,70),new=T) # Relabel.
#' }
#' @export fill.outside.border
fill.outside.border <-
function(col = 0, rat = 1)
{
	geopar <- getOption("geopar")
	gx <- geopar$limx
	gy <- geopar$limy
	gx <- mean(gx) + rat * (gx - mean(gx))
	gy <- mean(gy) + rat * (gy - mean(gy))
	dx <- gx[2] - gx[1]
	dy <- gy[2] - gy[1]
	x1 <- gx[1] - dx
	x2 <- gx[2] + dx
	y1 <- gy[1] - dy
	y2 <- gy[2] + dy
	b1 <- list(x = c(x1, x2, x2, x1, x1), y = c(gy[2], gy[2], y2, y2, gy[
		2]))
	b2 <- list(x = c(x1, x2, x2, x1, x1), y = c(gy[1], gy[1], y1, y1, gy[
		1]))
	b3 <- list(x = c(gx[2], x2, x2, gx[2], gx[2]), y = c(gy[1], gy[1],
		gy[2], gy[2], gy[1]))
	b4 <- list(x = c(gx[1], x1, x1, gx[1], gx[1]), y = c(gy[1], gy[1],
		gy[2], gy[2], gy[1]))
	oldpar <- selectedpar()
	par(geopar$gpar)
	polygon(b1, col = 0)
	polygon(b2, col = 0)
	polygon(b3, col = 0)
	polygon(b4, col = 0)
	par(oldpar)
	return(invisible())
}

