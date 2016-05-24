#' Aggregate data on square grid, pad with zeros.  Faster than tapply.
#' 
#' Function that smooths data so data.  In each given square where the number
#' of points exceeds certain minimum the program gives the mean, sum, variance
#' or median of the points within the square as well as the number of points
#' and indication if the minimum number of points was reached.  For each
#' datapoint latitude, longitude and z has to be given.  The program returns
#' the mean and mean longitude for points within a square, as well as the
#' number of points behind that point. The program also has the possibility to
#' remove outliers within each square. This program can be used to prepare data
#' for the program pointkriging if the data is so dense in certain places that
#' it disturbs the neighbourhood search or if there are +10000 datapoints.
#' 
#' It can also be used to find sums or averages of points within squares.
#' 
#' The program works quite fast even for very large datasets.
#' 
#' 
#' @param lat Latitude of datapoints. ( or x coordinates)
#' @param lon Longitude of datapoints. (or y coordinate)
#' @param z Values at datapoints.
#' @param grlat Latitude of defined grid.
#' @param grlon Longitude of defined grid.  If 0 grlat is a list with
#' components \$lat and \$lon.
#' @param fun Function (?).
#' @param fill If fill is T the center points of squares where the number of
#' points does not reach the minimum is set to zero.  Good for data like
#' trawlers report where it can be assumed that there is little fish where they
#' do not try.
#' @param reg List with two components \$lat and \$lon.  Points outside region
#' are not returned.  The list is typically output from the program
#' define.area.
#' @param minnumber Minimum number of points needed for the square to be valid.
#' Default value is 2.
#' @param wsp Workspace.  Not of interest.
#' @param wz Weighing for the z values.  If one is trying to find the mean
#' tons/hour in an area is it reasonable to weigh tons/hour by the catch.  If
#' the weight is zero within a square unweighed mean is used.
#' @param wlat Weighing to find the mean latitude and mean longitude.  It is
#' often reasonable to weigh the latitude and the longitude by the catch.
#' @param xy Projected? Default unprojected (lat/lon) values, \code{FALSE}.
#' @param rat Ratio for "rm.outliers".  If rat is for example 0.2 20\% of the
#' points on each side are moved.  i.e. the middle 60\% are kept.  Default
#' value is 0.6.
#' @param type Only value with description \code{rm.outliers}.
#' @section Value: <s-example> List with components \$lat,\$lon,\$z,\$n and
#' \$fill.  \$lat mean latitude of points within each square.  \$lon mean
#' longitude of points within each square \$z sum, mean, variance or median of
#' data \$n Number of points behind point.
#' 
#' If type = "rm.outliers" \$n this component stores the number of the square
#' corresponding to each point. \$fill If 0 there are datapoints behind the
#' points else it is an artificial point. (fill = T) </s-example>
#' @seealso \code{\link{variogram}}, \code{\link{variofit}},
#' \code{\link{pointkriging}}, \code{\link{grid}},
#' \code{\link{geocontour.fill}}.
#' @keywords <!--Put one or more s-keyword tags here-->
#' @export combine.rt
combine.rt <-
function(lat, lon, z, grlat, grlon = 0, fun, fill = F, reg = 0, minnumber = 2,
	wsp = 0, wz = 0, wlat = 0, xy = F, rat = 0.2, type)
{
	if(missing(fun) && !missing(type))
		fun <- type
	# for compatibility
	if(!missing(fun) && fun == "summa") fun <- "sum"
	# also for compatibility
	if(xy) {
		if(length(grlon) < 2) {
			grlon <- grlat$y
			grlat <- grlat$x
		}
	}
	else {
		if(length(grlon) < 2) {
			grlon <- grlat$lon
			grlat <- grlat$lat
		}
	}
	ndata <- length(lat)
	if(length(wz) != ndata)
		wz <- rep(1, ndata)
	if(length(wlat) != ndata)
		wlat <- rep(1, ndata)
	n <- length(grlon)
	m <- length(grlat)
	row <- cut(lat, grlat,labels=FALSE)  # R ver
	col <- cut(lon, grlon,labels=FALSE)  # R ver
	reitur <- (n - 1) * (row - 1) + col
	ind <- c(1:length(reitur))
	ind <- ind[!is.na(reitur)]
	lat <- lat[ind]
	lon <- lon[ind]
	z <- z[ind]
	wlat <- wlat[ind]
	wz <- wz[ind]
	reitur <- reitur[ind]
	maxrt <- (n - 1) * (m - 1)
	grdlat <- (grlat[1:(m - 1)] + grlat[2:m])/2
	grdlon <- (grlon[1:(n - 1)] + grlon[2:n])/2
	#	what to do 
	if(fun == "mean") option <- 1
	if(fun == "sum")
		option <- 2
	if(fun == "median")
		option <- 3
	if(fun == "variance")
		option <- 4
	if(fun == "rm.outliers")
		option <- 5
	if(fun == "keep.all")
		option <- 6
	#	Fill up matrix of data.  
	pts.in.reit <- c(matrix(0, round(ndata * 1.2), 1))
	npts.in.reit <- jrt <- indrt <- rep(0, maxrt + 1)
	nnewlat <- 0
	if(option == 5) {
		newlat <- newlon <- newz <- newn <- fylla <- rep(0, length(
			lat))
	}
	else if(option == 6) {
		newlat <- newlon <- newz <- newn <- fylla <- rep(0, (length(
			lat) + length(grlat) * length(grlon)) * 1.1)
	}
	else {
		newlat <- newlon <- newz <- newn <- fylla <- rep(0, maxrt +
			1)
	}
	if(wsp == 0)
		wsp <- ndata
	workspace <- order <- nr <- rep(0, wsp)
	outcome <- .C("combinert", PACKAGE = "geo", 
		as.double(lat),
		as.double(lon),
		as.double(z),
		as.integer(length(z)),
		as.integer(reitur),
		as.integer(pts.in.reit),
		as.integer(npts.in.reit),
		as.integer(indrt),
		as.integer(jrt),
		as.integer(maxrt),
		as.double(newlat),
		as.double(newlon),
		as.double(newz),
		as.integer(newn),
		as.integer(nnewlat),
		as.integer(option),
		as.integer(fill),
		as.integer(fylla),
		as.double(grdlat),
		as.double(grdlon),
		as.integer(n),
		as.integer(minnumber),
		as.double(workspace),
		as.integer(nr),
		as.integer(order),
		as.double(wlat),
		as.double(rat),
		as.double(wz))
	nnewlat <- outcome[[15]]
	newn <- outcome[[14]][1:nnewlat]
	newlat <- outcome[[11]][1:nnewlat]
	newlon <- outcome[[12]][1:nnewlat]
	newz <- outcome[[13]][1:nnewlat]
	fylla <- outcome[[18]][1:nnewlat]
	if(xy)
		projection <- "none"
	else projection <- "Mercator"
	if(length(reg) > 1) {
		inni <- inside(newlat, newlon, reg, option = 0, projection = 
			projection)
		ind <- c(1:length(inni))
		ind <- ind[inni == 1]
		newlat <- newlat[ind]
		newlon <- newlon[ind]
		newz <- newz[ind]
		newn <- newn[ind]
	}
	if(option == 5)
		fylla <- 0
	# not used
	if(xy) z <- list(x = newlat, y = newlon, z = newz, n = newn, fill = 
			fylla) else z <- list(lat = newlat, lon = newlon, z = 
			newz, n = newn, fill = fylla)
	z <- data.frame(z)
	attributes(z)$fun <- fun
	return(invisible(z))
}

