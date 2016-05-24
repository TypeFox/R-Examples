#' Shade border ?
#' 
#' Shade border ?
#' 
#' 
#' @param reg Region
#' @param lat Latitude ?
#' @param lon Longitude ?
#' @param col Color, not used ?
#' @param col.names Column names containing coordinates, default \code{lat} and
#' \code{lon}.
#' @return No value, addes shadedborder (with call to \code{\link{lines}}) to
#' contoured geoplot.
#' @note Needs elaboration. Color argument not used, but fixed values given to
#' 2 calls to \code{lines}.
#' @seealso Called by \code{\link{geocontour}}, calls \code{\link{Proj}}.
#' @keywords aplot
#' @export shadeborder
shadeborder <-
function(reg, lat, lon, col = 0, col.names = c("lon", "lat"))
{
	geopar <- getOption("geopar")
	ind <- c(1:length(reg[[col.names[2]]]))
	ind1 <- ind[is.na(reg[[col.names[2]]])]
	if(length(ind1) == 0 || ind1[1] != 1) {
		#external border does not begin with NA
		if(length(ind1) < 1) ind2 <- length(reg[[col.names[2]]]) else 
				ind2 <- ind1[1] - 1
		reg.lat <- reg[[col.names[2]]][1:ind2]
		reg.lon <- reg[[col.names[1]]][1:ind2]
		lonx <- c(min(lon), min(lon), max(lon), max(lon), min(lon),
			min(lon))
		latx <- c(mean(lat), min(lat), min(lat), max(lat), max(lat),
			mean(lat))
		ind2 <- ind[reg.lon == min(reg.lon)][1]
		ind3 <- ind[reg.lon == max(reg.lon)][1]
		ind6 <- ind[reg.lat == min(reg.lat)][1]
		ind7 <- ind[reg.lat == max(reg.lat)][1]
		i <- 0
		if(ind6 > ind2)
			i <- i + 1
		if(ind3 > ind6)
			i <- i + 1
		if(ind7 > ind3)
			i <- i + 1
		if(ind2 > ind7)
			i <- i + 1
		if(i > 1)
			ccw <- T
		else ccw <- F
		if(ccw) {
			#counterclockwise
			if(ind3 > ind2) {
				ind4 <- c(ind3:ind2)
				ind5 <- c(ind3:length(reg.lat), 1:ind2)
			}
			else {
				ind4 <- c(ind3:1, length(reg.lat):ind2)
				ind5 <- c(ind3:ind2)
			}
		}
		else {
			#clockwise
			if(ind3 > ind2) {
				ind4 <- c(ind3:length(reg.lat), 1:ind2)
				ind5 <- c(ind3:ind2)
			}
			else {
				ind4 <- c(ind3:ind2)
				ind5 <- c(ind3:1, length(reg.lat):ind2)
			}
		}
		mil <- min(min(lon), min(reg.lon) - 1)
		mal <- max(max(lon), max(reg.lon) + 1)
		rlat <- c(mean(lat), min(lat), min(lat), mean(lat))
		rlon <- c(mil, mil, mal, mal)
		rlon <- c(reg.lon[ind4], rlon)
		rlat <- c(reg.lat[ind4], rlat)
		rx <- Proj(rlat, rlon, geopar$scale, geopar$b0, geopar$b1,
			geopar$l1, geopar$projection, col.names = col.names)
		lines(rx, lwd = 2)
		#    polygon(rx$x, rx$y, border = F, col = col)
		rlat <- c(mean(lat), max(lat), max(lat), mean(lat))
		rlon <- c(mil, mil, mal, mal)
		rlon <- c(reg.lon[ind5], rlon)
		rlat <- c(reg.lat[ind5], rlat)
		rx <- Proj(rlat, rlon, geopar$scale, geopar$b0, geopar$b1,
			geopar$l1, geopar$projection, col.names = col.names)
		lines(rx, lwd = 2, col = 70)
		#    polygon(rx$x, rx$y, border = F, col = col)
		if(length(ind1) > 0) {
			if(geopar$projection == "none") {
				if(length(reg$x) - ind1[length(ind1)] < 3)
					return(invisible())
				reg$x <- reg$x[(ind1[1] + 1):length(reg$x)]
				reg$y <- reg$y[(ind1[1] + 1):length(reg$y)]
			}
			else {
				if(length(reg[[col.names[2]]]) - ind1[length(
					ind1)] < 3)
					return(invisible())
				reg[[col.names[2]]] <- reg[[col.names[2]]][
					(ind1[1] + 1):length(reg[[col.names[
					2]]])]
				reg[[col.names[1]]] <- reg[[col.names[1]]][
					(ind1[1] + 1):length(reg[[col.names[
					1]]])]
			}
			rx <- Proj(reg, scale = geopar$scale, b0 = geopar$
				b0, b1 = geopar$b1, l1 = geopar$l1, projection
				 = geopar$projection, col.names = col.names)
			lines(rx, lwd = 2, col = 2)
		}
	}
	else {
		rx <- Proj(reg, scale = geopar$scale, b0 = geopar$b0, b1 = 
			geopar$b1, l1 = geopar$l1, projection = geopar$
			projection, col.names)
		lines(rx, lwd = 2, col = 150)
	}
}

