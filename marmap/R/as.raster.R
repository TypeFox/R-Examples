as.raster <- function(bathy) {
	
	if (!is(bathy,"bathy")) stop("Objet is not of class bathy")

	lat <- as.numeric(colnames(bathy))
	lon <- as.numeric(rownames(bathy))

	r <- raster::raster(ncol = nrow(bathy), nrow = ncol(bathy), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
	raster::values(r) <- as.vector(bathy[,rev(1:ncol(bathy))])

	return(r)

}
