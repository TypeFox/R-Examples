#' Spherical (?) distance
#' 
#' Spherical (?) distance.
#' 
#' 
#' @param lat,lon Coordinate vectors of starting positions in lat/lon
#' @param lat1,lon1 Coordinate vectors of end positions in lat/lon
#' @return Distance in kilometers.
#' @note Why not use \code{arcdist}? Distance functins could be documented in
#' one file.
#' @seealso Called by \code{\link{pointkriging}} and \code{\link{variogram}}.
#' @keywords manip
#' @export pdist
pdist <-
function(lat, lon, lat1, lon1)
{
	rad <- 6367
	#radius of earth in km
	return(rad * acos(sin(lat) * sin(lat1) + cos(lat) * cos(lat1) * cos(
		lon - lon1)))
}

