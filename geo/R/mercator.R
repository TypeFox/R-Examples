#' Mercator projection
#' 
#' Mercator projection.
#' 
#' 
#' @param lat,lon Coordinates as latitude and longitude vectors
#' @param scale Scale of the output, "km" default, all other values imply
#' nautical miles.
#' @param b0 Latitude defining the projection.
#' @return List of components: \item{lat, lon }{Coordinates in latitude and
#' longitude} \item{x, y}{Input coordinates (projected)} \item{scale}{Scale}
#' \item{projection}{Projection (stated redunantly)} \item{b0, L0}{Defining lat
#' and a null value ???} is returned invisibly.
#' @note Needs elaboration, could/should (?) be documented with other
#' projection functions.
#' @seealso Called by \code{\link{geoarea}}, \code{\link{gridaxes}} and
#' \code{\link{Proj}}.
#' @keywords manip
#' @export mercator
mercator <-
function(lat, lon, scale = "km", b0 = 65)
{
	radius <- 6378.388
	m.p.km <- 1.852
	mult <- radius
	if(scale != "km")
		mult <- mult/m.p.km
	l1 <- (lon * pi)/180
	b1 <- (lat * pi)/180
	b0 <- (b0 * pi)/180
	x <- mult * cos(b0) * l1
	y <- mult * cos(b0) * log((1 + sin(b1))/cos(b1))
	return(invisible(list(lat = lat, lon = lon, x = x, y = y, scale = scale,
		projection = "mercator", b0 = b0, L0 = NULL)))
}

