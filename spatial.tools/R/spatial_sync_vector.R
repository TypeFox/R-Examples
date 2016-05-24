#' Matches a vector's projection to another vector or raster object's projection.
#' @param unsynced The vector to be projected.
#' @param reference A raster or vector object who's projection the unsynced will be matched to.
#' @param verbose Logical. Verbose logging?
#' @name spatial_sync_vector
#' @author Jonathan A. Greenberg (\email{spatial.tools@@estarcion.net})
#' @examples
#' tahoe_highrez_training_points_utm <- readOGR(
#' 	dsn=system.file("external", package="spatial.tools"),
#' 	layer="tahoe_highrez_training_points_utm")
#' print(projection(tahoe_highrez_training_points_utm))
#' tahoe_lidar_bareearth <- 
#' 	raster(system.file("external/tahoe_lidar_bareearth.tif", package="spatial.tools"))
#' print(projection(tahoe_lidar_bareearth))
#' tahoe_highrez_training_points_utm_synced <- 
#' 	spatial_sync_vector(tahoe_highrez_training_points_utm,tahoe_lidar_bareearth)
#' print(projection(tahoe_highrez_training_points_utm_synced))
#' @import rgdal
#' @export

spatial_sync_vector <- function(unsynced,reference,verbose=TRUE)
{
#	require("rgdal")
	new_projection=projection(reference)
	old_projection=projection(unsynced)

	if(new_projection!=old_projection)
	{
		synced_vector=spTransform(unsynced,CRS(new_projection))
	} else
	{
		synced_vector=unsynced
	}
	if(verbose){ message(projection(synced_vector)) } 
	return(synced_vector)
}