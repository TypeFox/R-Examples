#' Plot tows as line segments
#' 
#' The function gets 4 arguments i.e position of begininning and end of
#' segments.  If the first argument is given col.names gives the column names
#' in the data frame describing the position.
#' 
#' 
#' @param lat Vector of start position latitudes, dataframe of start positions 
#'  given in columns \code{lat} and \code{lon} or dataframe with start and end 
#'  positions given as columns \code{col.names}.
#' @param lon When given, vector of start position longitudes. 
#' @param lat1 When given, vector of end position latitudes or a data frame 
#'  of end positions given in columns \code{lat} and \code{lon}.
#' @param lon1 When given, vector of end position longitudes.
#' @param col Color of line segments.
#' @param col.names  Column names of start lat and lon and end lat and lon 
#'  when  \code{tows} are given in a single dataframe. Defaults to column
#'  names in the Hafro/MRI database table \code{fiskar.stodvar}.
#' @param \dots Additional arguments to \code{geolines}.
#' @seealso \code{\link{geolines}}
#' @keywords aplot
#' @export geotows
geotows <-
function(lat, lon, lat1, lon1, col=1,col.names = c("kastad.n.breidd", 
	"kastad.v.lengd", "hift.n.breidd", "hift.v.lengd"), ...)
{
	if(is.data.frame(lat) && missing(lat1)) {
		lat1 <- lat[, col.names[3]]
		lon1 <- lat[, col.names[4]]
		lon <- lat[, col.names[2]]
		lat <- lat[, col.names[1]]
	}
	if(is.data.frame(lat) && !missing(lat1)) {
		lon <- lat$lon
		lat <- lat$lat
		lon1 <- lat1$lon
		lat1 <- lat1$lat
	}
	lat <- matrix(lat, length(lat), 3)
	lat[, 2] <- lat1
	lat[, 3] <- NA
	lat <- c(t(lat))
	lon <- matrix(lon, length(lon), 3)
	lon[, 2] <- lon1
	lon[, 3] <- NA
	lon <- c(t(lon))
	geolines(lat, lon,col=col, ...)
	return(invisible())
}

