#' Fill points (thicken)
#' 
#' Fill points (thicken) for drawing continous lines in Lambert projection.
#' 
#' 
#' @param x,y Coordinates
#' @param nx Thickening factor
#' @param option Deals with NAs in the coordinates when not 1, the default
#' @return List of thickened values with components: \item{x, y}{of
#' coordinates}
#' @note Internal, needs elaboration.
#' @seealso The function is called by \code{\link{geopolygon}},
#' \code{\link{geolines}}, \code{\link{reitaplott}} and
#' \code{\link{gridaxes.Lambert}}
#' @keywords manip
#' @export fill.points
fill.points <-
function(x, y, nx, option = 1)
{
	n <- length(x)
	ny <- nx
	if(option != 1) {
		naind <- c(1:length(x))
		naind <- naind[is.na(x)]
	}
	dx <- (x[2:n] - x[1:(n - 1)])/(ny)
	dy <- (y[2:n] - y[1:(n - 1)])/(ny)
	x1 <- matrix(x[1:(n - 1)], n - 1, nx)
	y1 <- matrix(y[1:(n - 1)], n - 1, nx)
	ind <- c(0:(nx - 1))
	ind <- matrix(ind, n - 1, nx, byrow = T)
	dx <- matrix(dx, n - 1, nx)
	dy <- matrix(dy, n - 1, nx)
	x1 <- t(x1 + ind * dx)
	y1 <- t(y1 + ind * dy)
	ind <- c(1:length(y1))
	ind <- ind[is.na(y1) & row(y1) != 1]
	if(length(ind) != 0) {
		x1 <- x1[ - ind]
		y1 <- y1[ - ind]
	}
	if(is.na(x1[length(x1)])) {
		x1 <- c(x1, NA)
		y1 <- c(y1, NA)
	}
	ind <- c(1:length(x1))
	ind <- ind[is.na(x1)]
	if(length(ind) > 0) {
		ind <- matrix(ind,  , 2, byrow = T)
		if(option == 1) {
			ind <- ind[, 1]
			x1 <- x1[ - ind]
			y1 <- y1[ - ind]
		}
		else {
			ind <- ind[, 1]
			x1[ind] <- x[naind - 1]
			y1[ind] <- y[naind - 1]
		}
	}
	if(option != 1) {
		x1 <- c(x1, x[n])
		y1 <- c(y1, y[n])
	}
	return(list(x = x1, y = y1))
}

