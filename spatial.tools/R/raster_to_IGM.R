#' Creates an input geometry (IGM) file from a raster
#' 
#' @param x Raster* The raster containing the geographic information to be used as the basis for the IGM.
#' @author Jonathan A. Greenberg
#' @seealso \code{\link[raster]{xyFromCell}}, \code{\link[raster]{cellFromRow}}
#' 
#' @details This function produces a two-band brick where the pixel
#' values for the first band are the geographic x-coordinates, and the
#' pixel values for the second band are the geographic y-coordinates
#' of the corresponding pixels in the input Raster* file.
#' 
#' @import raster
#' @import foreach
#' @examples
#' tahoe_lidar_highesthit <- 
#' 	raster(system.file("external/tahoe_lidar_highesthit.tif", package="spatial.tools"))
#' tahoe_lidar_highesthit_igm <- raster_to_IGM(tahoe_lidar_highesthit)
#' plot(tahoe_lidar_highesthit_igm)
#' setMinMax(tahoe_lidar_highesthit_igm)
#' 
#' @export

raster_to_IGM <- function(x)
{
	igm_function <- function(x)
	{
		coords <- xyFromCell(x,cellFromRow(x,1:nrow(x)))
		dim(coords) <- c(dim(x)[2:1],2)
		return(coords)
	}
	
	igm <- rasterEngine(x=x,fun=igm_function,chunk_format="raster",debugmode=FALSE)
	return(igm)
}