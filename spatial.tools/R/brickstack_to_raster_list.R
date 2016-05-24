#' Converts a RasterBrick/RasterStack to a list of RasterLayers
#' @param x A RasterBrick or RasterStack.
#' @return A list of RasterLayers.
#' @author Jonathan A. Greenberg
#' @examples
#' # You can speed this up if a parallel backend is running, e.g.:
#' # sfQuickInit()
#' registerDoSEQ() # Just to avoid the warning from foreach.
#' tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' tahoe_highrez_list <- brickstack_to_raster_list(tahoe_highrez)
#' tahoe_highrez_list
#' # sfQuickStop()
#' @import foreach foreach
#' @import raster
#' @export

brickstack_to_raster_list <- function(x)
{
	layer <- NULL # To make CRAN happy
	nlayers_x <- nlayers(x)
	raster_list <- foreach(layer=1:nlayers_x,.packages=c("raster")) %dopar%
		raster(x,layer=layer)
	
	return(raster_list)
	
}