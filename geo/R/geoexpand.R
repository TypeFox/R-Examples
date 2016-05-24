#' expands a grid to a dataframe.
#' 
#' Expands a grid to a dataframe simular to expand.grid.
#' 
#' 
#' @param grid The grid to be expanded, a list containing vectors \$lat and
#' \$lon.
#' @return A dataframe of the grid.
#' @section Side Effects: None.
#' @seealso \code{\link{expand.grid}}.
#' @examples
#' 
#'        grd <- list(lat=seq(63,67,length=30),lon=seq(-28,-10,length=50))
#'        # a list with length(lat) = 30 and length(lon) = 50
#' 
#'        grd1<- geoexpand(grd)
#'        # a 30 by 50 dataframe made.
#' 
#'        # See also examples in geocontour.
#' 
#' @export geoexpand
geoexpand <-
function(grid)
{
	if(is.null(grid$lat)) {
		ny <- length(grid$y)
		nx <- length(grid$x)
		y <- t(matrix(grid$y, ny, nx))
		x <- matrix(grid$x, nx, ny)
		return(data.frame(y = c(y), x = c(x)))
	}
	else {
		nlat <- length(grid$lat)
		nlon <- length(grid$lon)
		lat <- t(matrix(grid$lat, nlat, nlon))
		lon <- matrix(grid$lon, nlon, nlat)
		return(data.frame(lat = c(lat), lon = c(lon)))
	}
}

