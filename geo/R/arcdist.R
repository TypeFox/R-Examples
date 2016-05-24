#' Geographical distance computations
#' 
#' Computes distances between lat/lon data points.
#' 
#' 
#' @param lat Latitude of first coordinate or list with lat, lon of first
#' coordinate.
#' @param lon Longitude of first coordinate or list with lat, lon of second
#' coordinate.
#' @param lat1,lon1 If lat and lon are vectors of lat,lon positions, then lat1
#' and lon1 must be given as the second set of positions.
#' @param scale \code{nmi} returns value in nautical miles, any other value in
#' kilometers
#' @return A single vector of distances between pairs of points is returned
#' @seealso \code{\link{geoplot}}, \code{\link{geotext}}, \code{\link{selpos}},
#' @examples
#' 
#'   pos1 <- list(lat = c(65, 66), lon = c(-19, -20))
#'   pos2 <- list(lat = c(64, 65), lon = c(-19, -20))
#'   dists <- arcdist(pos1, pos2)         # pos1 and pos2 are lists of coordinates.
#'   lat <- c(65, 66)
#'   lon <- c(-19, -20)
#'   lat1 <- c(64, 65)
#'   lon1 <- c(-19, -20)
#'   dists <- arcdist(lat, lon, lat1, lon1) # Input in vector format.
#' 
#' @export arcdist
arcdist <-
function(lat, lon, lat1 = NULL, lon1 = NULL, scale = "nmi")
{
	if(is.null(lat1)) {
		lat1 <- lon$lat
		lon1 <- lon$lon
		lon <- lat$lon
		lat <- lat$lat
	}
	if(scale == "nmi")
		miles <- 1.852
	else miles <- 1
	rad <- 6367
	#radius of earth in km
	mult1 <- (rad/miles)
	mult2 <- pi/180
	return(mult1 * acos(sin(mult2 * lat) * sin(mult2 * lat1) + cos(mult2 *
		lat) * cos(mult2 * lat1) * cos(mult2 * lon - mult2 * lon1)))
}

