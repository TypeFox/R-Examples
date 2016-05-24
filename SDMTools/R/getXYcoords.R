#' Computes the X and Y Coordinates of the Pixels of a Raster Map
#' 
#' \code{getXYcoords} computes the geographical coordinates of the rows and
#' columns of pixels of a raster map of class \code{asc}. Code & helpfile were
#' modified from adehabitat package.
#' 
#' 
#' @param w an object of class \code{asc}.
#' @return Returns a list with two components: \item{x}{the x coordinates of
#' the columns of pixels of the map} \item{y}{the y coordinates of the rows of
#' pixels of the map}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' 
#' 
#' tasc = as.asc(matrix(rep(x=1:10, times=1000),nr=100)); print(tasc)
#' getXYcoords(tasc)
#' 
#' 
#' @export 
"getXYcoords" <- function(w)
{
	#check if raster from sp or raster package and convert if necessary
	if (any(class(w) %in% 'RasterLayer')) w = asc.from.raster(w)
	if (any(class(w) == 'SpatialGridDataFrame')) w = asc.from.sp(w)
    if (!inherits(w, "asc")) stop("must be of class asc")

    # Gets the attributes
    cs<-attr(w, "cellsize")
    xll<-attr(w, "xll")
    yll<-attr(w, "yll")

    ## Computation of the number of rows and columns of the matrix
    nr<-nrow(w)
    nc<-ncol(w)

    ## The results
    x<-xll+c(0:(nr-1))*cs
    y<-yll+c(0:(nc-1))*cs
    return(list(x=x, y=y))
}

