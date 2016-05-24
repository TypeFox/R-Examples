#' Zoom into plots.
#' 
#' Zoom into current geoplot, redraws current geoplot with area defined by
#' user.  The current plot can be restored by using geodezoom.
#' 
#' 
#' @return none.
#' @section Side Effects: The last call to geoplot is recalled with different
#' borders for x and y.
#' @seealso \code{\link{geoplot}}, \code{\link{geodezoom}}.
#' @examples
#' 
#'     geoplot()
#'     geozoom()
#'     # Click with mouse as when placing legend i.e. first place the
#'     # mouse where the upper left corner is supposed to be and click
#'     # once, then move to where the lower right corner is supposed to
#'     # be and also press once.
#'     geodezoom()
#'     # Return to previous plot, here geoplot().
#' 
#' @export geozoom
geozoom <-
function()
{
	geopar <- getOption("geopar")
	if(as.character(geopar$command[length(geopar$command)]) != "123")
		com <- c(geopar$command, zoom = 123)
	else com <- geopar$command
	eval(com)
}

