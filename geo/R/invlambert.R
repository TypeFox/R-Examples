#' Inverse Lambert projection
#' 
#' Inverse Lambert projection.
#' 
#' 
#' @param x,y The input data to be inversely projected as two vectors
#' @param lat0 First latitude defining of the projection
#' @param lon0 Longitude defining the projection
#' @param lat1 Second latitude defining the projection
#' @param scale Scale of the input data, "km" default, redundant ??
#' @param old Old method, seldom used???
#' @return List of components: \item{lat, lon }{Coordinates in latitude and
#' longitude} \item{x, y}{Input coordinates (projected)} \item{scale}{Scale}
#' \item{projection}{Projection (stated redunantly)} \item{lat0, lon0,
#' lat1}{Defining lats and lon} is returned invisibly.
#' @note Needs elaboration and perhaps documenting with lambert in the same
#' doc-file.
#' @seealso Called by \code{\link{invProj}}.
#' @keywords manip
#' @export invlambert
invlambert <-
function(x, y, lat0, lon0, lat1, scale = "km", old = F)
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
	p <- sign(n) * sqrt(x^2 + (p0 - y)^2)
	theta <- atan(x/(p0 - y))
	t <- (p/(a * F1))^(1/n)
	lon <- theta/n + lon0
	lat <- pi/2 - 2 * atan(t)
	for(i in 1:2) {
		# very rapid convergence in all cases.  
		lat <- pi/2 - 2 * atan(t * ((1 - e * sin(lat))/(1 + e * sin(
			lat)))^(e/2))
	}
	return(invisible(list(lat = (lat * 180)/pi, lon = (lon * 180)/pi, x = x,
		y = y, scale = scale, projection = "lambert", lat0 = (lat0 *
		180)/pi, lon0 = (lon0 * 180)/pi, lat1 = lat11)))
}

