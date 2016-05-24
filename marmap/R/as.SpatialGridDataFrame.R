as.SpatialGridDataFrame <- function(bathy) {
	
	out <- as(as.raster(bathy), "SpatialGridDataFrame")
	return(out)
}


