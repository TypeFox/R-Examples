#' Convert latitude and longitude
#' 
#' Convert between different representations of latitude and longitude, namely
#' degrees-minutes-decimal minutes and decimal degrees.
#' 
#' 
#' @param data Dataframe with coordinates in two columns
#' @param inverse Which conversion should be undertaken, default from
#' degrees-minutes-decimal minutes (DDMMmm) to decimal degrees (DD.dd)
#' @param col.names Colnames in \code{data} with coordinates to convert,
#' default \code{lat, lon}
#' @return Returns \code{data} with converted values in the coordinate columns.
#' @note Functions calling \code{geoconvert} do so in a branch that is probably
#' rarely used. Implement conversion from other representations of lat and lon
#' in future?
#' @seealso Called by a number of functions, i.e. \code{\link{d2mr}},
#' \code{\link{mr2d}}, \code{\link{geogrid}}, \code{\link{geoidentify}},
#' \code{\link{geolines}}, \code{\link{geopoints}}, \code{\link{geopolygon}}
#' and \code{\link{geotext}}, perhaps unncessarily in some. Wraps around
#' functions \code{\link{geoconvert.1}} and \code{\link{geoconvert.2}},
#' depending on which conversion to undertake.
#' @keywords manip
#' @export geoconvert
geoconvert <-
function(data, inverse = F, col.names = c("lat", "lon"))
{
	if(!inverse) {
		if(is.data.frame(data)) {
			if(any(is.na(match(col.names, names(data))))) {
				cat(paste("Columns", colnames, "do not exist"))
				return(invisible())
			}
			data[, col.names[1]] <- geoconvert.1(data[, col.names[
				1]])
			data[, col.names[2]] <- geoconvert.1(data[, col.names[
				2]])
		}
		else data <- geoconvert.1(data)
	}
	else {
		# Convert to write out. 
		if(is.data.frame(data)) {
			if(any(is.na(match(col.names, names(data))))) {
				cat(paste("Columns", colnames, "do not exist"))
				return(invisible())
			}
			data[, col.names[1]] <- geoconvert.2(data[, col.names[
				1]])
			data[, col.names[2]] <- geoconvert.2(data[, col.names[
				2]])
		}
		else data <- geoconvert.2(data)
	}
	return(data)
}

