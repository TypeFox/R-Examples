#' Lambert projection
#' 
#' Lambert projection.
#' 
#' 
#' @param lat,lon Coordinates as latitude and longitude vectors
#' @param lat0 First latitude defining of the projection
#' @param lon0 Longitude defining the projection
#' @param lat1 Second latitude defining the projection
#' @param scale Scale, default "km", redundant ??
#' @param old Old method ?
#' @return List of components: \item{lat, lon }{Coordinates in latitude and
#' longitude} \item{x, y}{Input coordinates (projected)} \item{scale}{Scale}
#' \item{projection}{Projection (stated redunantly)} \item{lat0, lon0,
#' lat1}{Defining lats and lon} is returned invisibly.
#' @note Needs elaboration, might be merged with docs for other proj-functions.
#' @seealso Called by \code{\link{geoarea}}, \code{\link{orthproj}} and
#' \code{\link{Proj}}.
#' @keywords manip
#' @export lambert
lambert <-
function(lat, lon, lat0, lon0, lat1, scale = "km", old = F)
{
	a <- 6378.388
	# radius at equator
	e <- sqrt(2/297 - (1/297)^2)
	# eccensitret.
	lat11 <- lat1
	# temporary storage
	# 	change to radians  
	lat1 <- (lat1 * pi)/180
	lat0 <- (lat0 * pi)/180
	lon0 <- (lon0 * pi)/180
	lat <- (lat * pi)/180
	lon <- (lon * pi)/180
	#	one or two touching points.  
	if(length(lat1) == 2) {
		lat2 <- lat1[2]
		lat1 <- lat1[1]
		np <- 2
	}
	else np <- 1
	m1 <- cos(lat1)/sqrt(1 - e * e * (sin(lat1))^2)
	if(old) {
		t1 <- tan(pi/4 - 1/2 * atan((1 - e * e) * tan(lat1)))
		t0 <- tan(pi/4 - 1/2 * atan((1 - e * e) * tan(lat0)))
	}
	else {
		t1 <- tan(pi/4 - lat1/2)/((1 - e * sin(lat1))/(1 + e * sin(
			lat1)))^(e/2)
		t0 <- tan(pi/4 - lat0/2)/((1 - e * sin(lat0))/(1 + e * sin(
			lat0)))^(e/2)
	}
	# one tangent.   
	if(np == 1) n <- sin(lat1) else {
		m2 <- cos(lat2)/(1 - e * e * (sin(lat2))^2)
		if(old)
			t2 <- tan(pi/4 - 1/2 * atan((1 - e * e) * tan(lat2)))
		else t2 <- tan(pi/4 - lat2/2)/((1 - e * sin(lat2))/(1 + e *
				sin(lat2)))^(e/2)
		n <- (log(m1) - log(m2))/(log(t1) - log(t2))
	}
	F1 <- m1/(n * t1^n)
	p0 <- a * F1 * t0^n
	if(old)
		t <- tan(pi/4 - 1/2 * atan((1 - e * e) * tan(lat)))
	else t <- tan(pi/4 - lat/2)/((1 - e * sin(lat))/(1 + e * sin(lat)))^
			(e/2)
	p <- a * F1 * t^n
	theta <- n * (lon - lon0)
	x <- p * sin(theta)
	y <- p0 - p * cos(theta)
	return(invisible(list(lat = (lat * 180)/pi, lon = (lon * 180)/pi, x = x,
		y = y, scale = scale, projection = "Lambert", lat0 = (lat0 *
		180)/pi, lon0 = (lon0 * 180)/pi, lat1 = lat11)))
}

