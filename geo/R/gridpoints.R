#' Produce gridpoints over an area
#' 
#' Produce gridpoints over an area.
#' 
#' 
#' @param border Border of the area
#' @param dx Resolution in each direction (?)
#' @param grpkt ???
#' @param nx Number of gridpoints in each direction (?)
#' @param n Total number of gridpoints (?)
#' @return List with components: \item{xgr}{List of gridpoints in components
#' \code{lat, lon} or \code{x,y} depending on the projection.} \item{xgra}{List
#' of those gridpoints within the area given in \code{border}}
#' @note Needs further elaboration, check use with \code{find = TRUE} in
#' \code{setgrid}.
#' @seealso Called by \code{\link{setgrid}}.
#' @keywords manip
#' @export gridpoints
gridpoints <-
function(border, dx, grpkt, nx, n)
{
	geopar <- getOption("geopar")
	if(length(grpkt) == 1) {
		# gridpoints not given.  
		if(geopar$projection == "none") {
			xmin <- min(border$x)
			xmax <- max(border$x)
			ymin <- min(border$y)
			ymax <- max(border$y)
		}
		else {
			xmin <- min(border$lon)
			xmax <- max(border$lon)
			ymin <- min(border$lat)
			ymax <- max(border$lat)
			meanlat <- (mean(border$lat) * pi)/180
		}
		if(dx[1] == 0) {
			if(nx[1] == 0) {
				n <- sqrt(n)
				if(geopar$projection == "none")
					k <- (xmax - xmin)/(ymax - ymin)
				else k <- ((xmax - xmin) * cos(meanlat))/(
						ymax - ymin)
				nx[1] <- round(n * sqrt(k))
				nx[2] <- round(n/sqrt(k))
			}
			dx[1] <- (xmax - xmin)/nx[1]
			dx[2] <- (ymax - ymin)/nx[2]
		}
		else {
			tmp <- dx[1]
			dx[1] <- dx[2]
			dx[2] <- tmp
			#exchange lat lon.  
			nx[1] <- trunc((xmax - xmin)/dx[1])
			nx[2] <- trunc((ymax - ymin)/dx[2])
		}
		xgr <- (xmin - dx[1]) + c(1:(nx[1] + 2)) * dx[1]
		ygr <- (ymin - dx[2]) + c(1:(nx[2] + 2)) * dx[2]
	}
	else if(geopar$projection == "none") {
		xgr <- grpkt$x
		ygr <- grpkt$y
	}
	else {
		xgr <- grpkt$lon
		ygr <- grpkt$lat
	}
	lx <- length(xgr)
	ly <- length(ygr)
	xgra <- c(matrix(xgr, lx, ly))
	ygra <- c(t(matrix(ygr, ly, lx)))
	if(geopar$projection == "none") {
		xgra <- list(x = xgra, y = ygra)
		xgr <- list(x = xgr, y = ygr)
	}
	else {
		xgra <- list(lon = xgra, lat = ygra)
		xgr <- list(lon = xgr, lat = ygr)
	}
	return(list(xgr = xgr, xgra = xgra))
}

