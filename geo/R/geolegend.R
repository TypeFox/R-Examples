#' Put a legend on a plot in the geo series.
#' 
#' Adds a legend to current plot.  The location can be specified with \$lat and
#' \$lon.  Allows all the same parameters as legend.
#' 
#' See legend.
#' 
#' @param pos The position of the text, should include \$lat and \$lon.  If
#' \$lat and \$lon are llength 1, they determine the top left corner of the
#' rectangle; if theey are length 2 vectors, the give opposite corners of the
#' rectangular area.  A list containing x and y values may be supplied.
#' @param legend Vector of character strings to be associated with plot.
#' @param \dots The function allows any optional argument to the legend
#' function to be taken.
#' @return None.
#' @section Side Effects: Draws a box at specified coordinates and puts inside
#' (if possible) examples of lin, points, marks and/or shading, each identified
#' with a user-specified text string.
#' @seealso \code{\link{legend}}.
#' @examples
#' 
#'            # The function is currently defined as
#'        function(pos, legend, ...)
#'        {
#'                       oldpar <- par()
#'                       par(geopar$gpar)
#'                       on.exit(par(oldpar))
#'                       xx <- Proj(pos$lat, pos$lon)
#'                       legend(xx$x, xx$y, legend = legend, ...)
#'        }
#' 
#' @export geolegend
geolegend <-
function(pos, legend, ...)
{
	geopar <- getOption("geopar")
	oldpar <- selectedpar()
	par(geopar$gpar)
	on.exit(par(oldpar))
	xx <- Proj(pos$lat, pos$lon)
	legend(xx$x, xx$y, legend = legend, ...)
}

