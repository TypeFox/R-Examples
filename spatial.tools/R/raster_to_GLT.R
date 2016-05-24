#' Creates a geometry lookup (GLT) file from a raster
#' 
#' @param x Raster* The raster containing the geographic information to be used as the basis for the GLT.
#' @author Jonathan A. Greenberg
#' @seealso \code{\link[raster]{rowColFromCell}}, \code{\link[raster]{cellFromRow}}
#' 
#' @details This function produces a two-band brick where the pixel
#' values for the first band are the column numbers, and the
#' pixel values for the second band are the row numbers
#' of the corresponding pixels in the input Raster* file.
#' 
#' @import raster
#' @import foreach
#' @examples
#' tahoe_lidar_highesthit <- 
#' 	raster(system.file("external/tahoe_lidar_highesthit.tif", package="spatial.tools"))
#' tahoe_lidar_highesthit_glt <- raster_to_GLT(tahoe_lidar_highesthit)
#' plot(tahoe_lidar_highesthit_glt)
#' setMinMax(tahoe_lidar_highesthit_glt)
#' 
#' @export

raster_to_GLT <- function(x)
{
	glt_function <- function(x)
	{
		coords <- rowColFromCell(x,cellFromRow(x,1:nrow(x)))
		coords <- cbind(coords[,2],coords[,1])
		dim(coords) <- c(dim(x)[2:1],2)
		return(coords)
	}
	
	glt <- rasterEngine(x=x,fun=glt_function,chunk_format="raster",debugmode=FALSE)
	return(glt)
}