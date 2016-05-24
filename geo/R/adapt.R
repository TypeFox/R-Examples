#' Adapts geographical positions
#' 
#' Adapts geographical positions for further geo-use.
#' 
#' 
#' @param reg.lat Latitude or y-coordinate
#' @param reg.lon Longitude or x-coordinate
#' @param projection Projection, default "Mercator", "none" denotes x/y
#' coordinates
#' @return Returns a list of either: \item{x, y}{x- and y-positions} or:
#' \item{lat, lon}{Latitude and longitude} and: \item{lxv}{index of
#' uninterupted/contiguous positions}
#' @note Needs further elaboration, this function is called by
#' \code{geoarea.old}, \code{geoinside}, \code{inside} and \code{pointkriging}.
#' @keywords manip
#' @export adapt
adapt <-
function(reg.lat, reg.lon, projection = "Mercator")
{
	ind <- c(1:length(reg.lat))
	nholes <- length(reg.lat[is.na(reg.lat)])
	lxv <- c(0:(nholes + 1))
	if(nholes != 0) {
		#               remove NA,s and points given twice
		ind1 <- ind[is.na(reg.lat)]
		ind2 <- c(ind1 - 1, ind1, length(reg.lat))
		lon <- reg.lon[ - ind2]
		lat <- reg.lat[ - ind2]
		for(i in 2:(nholes + 1)) {
			lxv[i] <- ind1[i - 1] - 2 * (i - 1)
		}
	}
	else {
		ind <- (1:length(reg.lon) - 1)
		lon <- reg.lon[ind]
		lat <- reg.lat[ind]
	}
	lxv[nholes + 2] <- length(lon)
	#x,y coordinates.   
	if(projection == "none") return(list(x = lon, y = lat, lxv = lxv))
		 else return(list(lat = lat, lon = lon, lxv = lxv))
}

