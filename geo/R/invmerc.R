#' Inverse Mercator Projection
#' 
#' Inverse Mercator Projection.
#' 
#' 
#' @param x,y The input data to be inversely projected as two vectors
#' @param scale Scale of the input data, "km" default, all other values imply
#' nautical miles.
#' @param b0 Latitude defining the projection.
#' @return List of components: \item{lat, lon }{Coordinates in latitude and
#' longitude} \item{x, y}{Input coordinates (projected)} \item{scale}{Scale}
#' \item{projection}{Projection (stated redunantly)} \item{b0, L0}{Defining lat
#' and a null value ???} is returned invisibly.
#' @note Needs elaboration and perhaps documenting with mercator in the same
#' doc-file.
#' @seealso Called by \code{\link{invProj}}.
#' @keywords manip
#' @export invmerc
invmerc <-
function(x, y, scale = "km", b0 = 65)
{
	radius <- 6378.388
	m.p.km <- 1.852
	mult <- radius
	if(scale != "km")
		mult <- mult/m.p.km
	b0 <- (b0 * pi)/180
	lon <- (x/(mult * cos(b0)) * 180)/pi
	# Have to find latitude by iteration.
	c1 <- exp(y/(mult * cos(b0)))
	lat1 <- c(1:length(y))
	lat1[1:length(y)] <- b0
	# initial guess
	lat <- c(1:length(y))
	ind <- c(1:length(y))
	#index.
        ind <- ind[!is.na(c1)]
	# NA dont work in sum.
	while(sum(abs(lat1[ind] - lat[ind]))/sum(abs(lat[ind])) > 1e-07) {
		lat <- lat1
		lat1 <- lat - ((1 + sin(lat))/cos(lat) - c1)/((1 + sin(lat))/
			(cos(lat)^2))
	}
	lat <- lat1
	lat <- (lat * 180)/pi
	return(invisible(list(lat = lat, lon = lon, x = x, y = y, scale = scale,
		projection = "mercator", b0 = b0, L0 = NULL)))
}

