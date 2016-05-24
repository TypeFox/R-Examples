#' Convert from decimal degrees
#' 
#' Convert from decimal degrees to degrees, minutes and fractional minutes
#' representation (DDMMmm) of lat or lon.
#' 
#' 
#' @param lat Vector of latitude or longitudes
#' @return Returns a vector of six digit values with degrees, minutes and
#' fractions of minutes, with two decimal values, concatenated.
#' @seealso Called by \code{\link{geoconvert}}, when \code{inverse = TRUE}.
#' @keywords manip
#' @export geoconvert.2
geoconvert.2 <-
function(lat)
{
	i <- sign(lat)
	lat <- abs(lat)
	p1 <- floor(lat)
	p2 <- floor((lat - p1) * 60)
	p3 <- round((lat - p1 - p2/60) * 100 * 60)
	return(i * (p1 * 10000 + p2 * 100 + p3))
}

