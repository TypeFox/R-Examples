#' Defines regions.
#' 
#' locates points pointed out by users and defines regions.
#' 
#' Draws the regions on the plot.
#' 
#' @param nholes The number of holes in data, number of regions - 1.
#' @return A list of the points pointed out by the user with NA's between
#' regions.
#' @section Side Effects: Draws the regions on the plot.
#' @seealso \code{\link{geolocator}}.
#' @examples
#' 
#'  ##   Push left mouse button to mark point, push middle button to 
#'  ##   mark the end of a region.
#' 
#' @export geodefine
geodefine <-
function(nholes = 0)
{
	geopar <- getOption("geopar")
	oldpar <- selectedpar()
	par(geopar$gpar)
	on.exit(par(oldpar))
	border <- giveborder(nholes = nholes)
	reg <- border$reg
	par(oldpar)
	if(geopar$projection == "none")
		reg <- list(x = reg$x, y = reg$y)
	else reg <- list(lat = reg$lat, lon = reg$lon)
	reg <- data.frame(reg)
	return(reg)
}

