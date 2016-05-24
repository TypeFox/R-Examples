#' Calculates the distance between each pair of datapoints.
#' 
#' The program calculates the distance between each pair of datapoints.  The
#' distances are grouped in groups such that even number of pairs is in each
#' group.  Then the estimated variogram for each group is calculated either by
#' taking the mean or by a method from Cressie and Hawkins (1980).  The latter
#' method does in essence take the sum of the values^0.25. Only pair of points
#' with distance less than certain distance are used. Zero - Zero pairs are not
#' used if zzp is F.
#' 
#' 
#' @param lat Latitude of datapoints.  If lon=0 lat is a list with components
#' \$lat and \$lon.
#' @param z Values at datapoints.
#' @param lon longitude of datapoints.
#' @param nbins Number of distance intervals used.
#' @param maxdist maximum distance of interest.  Default value is range of
#' data*0.7.
#' @param Hawk If true the method from Cressie and Hawkins (1980) is used, else
#' the mean.  Default value is T
#' @param throwout If T datapoints with value zero are not used at all.
#' Default value is F which means that they are used.  In all cases zero-zero
#' pairs are not used when the variogram is estimated.
#' @param scale "km" or "nmi". Default is "km".
#' @param evennumber If T distance classes are chosen so approximately the same
#' number of points is in each distance class.  Else even distance increments
#' are used.
#' @param zzp If true zero-zero pairs are used else not.  Default value is F.
#' @param minnumber Distance intervals with minnumber or less pairs are not
#' included.  Default is zero.
#' @param col.names if lat is a dataframe col.names should contain the names of
#' the vectors containing the x and y coordinates, default is c("lat","lon")
#' @section Value: A list with the following components: <s-example>
#' 
#' number: Number of pair in each distance class.  dist: Mean distance in each
#' distance class.  vario: variogram for each distance class. </s-example> The
#' list is suitable for the program variofit.  The variogram can also be
#' plotted by plvar(vagram,fit=F)
#' @seealso \code{\link{variofit}}, \code{\link{pointkriging}},
#' \code{\link{plvar}}.
#' @export variogram
variogram <-
function(lat, lon = 0, z, nbins = 100, maxdist = 0, Hawk = T, throwout = F,
	scale = "km", evennumber = T, zzp = F, minnumber = 0, col.names = c(
	"lat", "lon"))
{
	if(is.data.frame(lat)) {
		lon <- lat[[col.names[2]]]
		lat <- lat[[col.names[1]]]
	}
	eps <- 1e-06
	rad <- 6378.388
	# Radius of earth in km.
	if(col.names[1] == "lat" && col.names[2] == "lon") {
		xy <- F
		if(length(lon) < 2) {
			lon <- lat$lon
			lat <- lat$lat
		}
		#list  
		if(scale == "nmi") rad <- rad/1.852
		# distances in miles.  
		lon <- (lon * pi)/180
		# change from degrees to radians
		lat <- (lat * pi)/180
	}
	else xy <- T
	if(throwout) {
		# throw out zero points.  
		lat <- lat[abs(z) > eps]
		lon <- lon[abs(z) > eps]
		z1 <- z[abs(z) > eps]
		z <- z1
	}
	variance <- var(z)
	count <- length(lon)
	# measurements
	if(maxdist == 0) {
		rlat <- range(lat)
		rlon <- range(lon)
		if(xy)
			maxdist <- pdistx(rlat[2], rlon[2], rlat[1], rlon[
				1]) * 0.7
		else maxdist <- pdist(rlat[2], rlon[2], rlat[1], rlon[1]) *
				0.7
		if(scale == "nmi")
			maxdist <- maxdist/1.852
	}
	varioa <- dista <- numbera <- rep(0, nbins)
	if(evennumber)
		nbins <- nbins * 10
	ddist <- maxdist/nbins
	dist <- vario <- number <- rep(0, nbins)
	dist <- .C("variogram", PACKAGE = "geo", 
		as.double(lat),
		as.double(lon),
		as.double(z),
		as.double(ddist),
		as.integer(number),
		as.double(dist),
		as.double(vario),
		as.integer(count),
		as.integer(nbins),
		as.integer(Hawk),
		as.integer(evennumber),
		as.double(varioa),
		as.double(dista),
		as.integer(numbera),
		as.integer(zzp),
		as.integer(xy))
	nbins <- dist[[9]]
	vario <- dist[[7]]
	xh <- dist[[6]]
	number <- dist[[5]]
	if(evennumber) {
		vario <- vario[1:nbins]
		xh <- xh[1:nbins]
		number <- number[1:nbins]
	}
	dist <- xh
	ind <- c(1:length(number))
	ind <- ind[number < minnumber + 1]
	if(length(ind) == 0)
		return(list(vario = vario, dist = dist, number = number, 
			variance = variance))
	else return(list(vario = vario[ - ind], dist = dist[ - ind], number = 
			number[ - ind], variance = variance))
}

