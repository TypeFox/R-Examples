#' Tests if an input is a RasterLayer, RasterBrick, or a RasterStack.
#' 
#' @param x an R Object.
#' @return A logical vector.
#' @author Jonathan A. Greenberg
#' @examples
#' tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' is.Raster(tahoe_highrez)
#' tahoe_lidar_bareearth <- 
#' 	raster(system.file("external/tahoe_lidar_bareearth.tif", package="spatial.tools"))
#' is.Raster(tahoe_lidar_bareearth)
#' is.Raster("character")
#' @import raster
#' @export

is.Raster <- function(x)
{
	return((class(x)=="RasterLayer" || class(x)=="RasterBrick" || class(x)=="RasterStack"))
}