#' Spatial Join of Points with Raster Grids - replace data
#' 
#' \code{put.data} replaces data in raster object of class 'asc' (this and
#' adehabitat package) at specified locations.\cr \cr \bold{Note:} there is no
#' interpolation done here. The values given replace the values of the raster
#' cell the point falls into.
#' 
#' Implements a faster version of 'join.asc' from the adehabitat package. \cr
#' \cr \bold{NOTE:} this assumes all locations are within the extent of the
#' raster map. Values outside the extent will be given a value of NA.
#' 
#' @param pts a three-column data frame or matrix with the x and y coordinates
#' of the locations of interest and the third column being the z values to put
#' in the ascii grid file.
#' @param x a raster matrix of class 'asc' (this and the adehabitat package)
#' @return Returns a raster matrix of class 'asc' equal in size to input 'x'.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' 
#' 
#' #create a simple object of class 'asc'
#' tasc = as.asc(matrix(1:50,nr=50,nc=50)); print(tasc)
#' \dontrun{image(tasc)}
#' 
#' #define some point locations
#' points = data.frame(x=runif(25,1,50),y=runif(25,1,50),z=50)
#' 
#' #put the new data
#' tasc = put.data(points,tasc)
#' 
#' #show the data
#' \dontrun{image(tasc)}
#' 
#' 
#' @export put.data
put.data <-
function(pts, x) {
	if (class(x) != 'asc') stop('matrix must be of class "asc"') #check to ensure x is of class asc
	pts = as.matrix(pts); if (dim(pts)[2]!=3) stop('input pts must be 3 columns of data')
	xy <- getXYcoords(x)
	xy$x <- xy$x + attr(x, "cellsize")/2
	xy$x <- c(xy$x[1] - attr(x, "cellsize"),xy$x)
	xy$y <- xy$y + attr(x, "cellsize")/2
	xy$y <- c(xy$y[1] - attr(x, "cellsize"),xy$y)
	xf <- as.numeric(cut(pts[, 1], xy$x))
	yf <- as.numeric(cut(pts[, 2], xy$y))
	x[cbind(xf,yf)] = pts[,3]
	return(x)
}
