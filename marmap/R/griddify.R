griddify <- function(xyz, nlon, nlat) {
    
	
	if (ncol(xyz) != 3) stop("xyz must be a 3-column matrix or data.frame of longitudes, latitudes and depths/altitudes")
	
    colnames(xyz) <- c("lon", "lat", "depth")
    
    r <- raster::raster(xmn = min(xyz$lon, na.rm = TRUE),
                        xmx = max(xyz$lon, na.rm = TRUE),
                        ymn = min(xyz$lat, na.rm = TRUE),
                        ymx = max(xyz$lat, na.rm = TRUE), ncol = nlon, nrow = nlat)
    
    x <- raster::rasterize(xyz[, c("lon", "lat")], r, xyz$depth, fun = mean)
    
    s <- raster::raster(nrow = nlat, ncol = nlon)
    raster::extent(s) <- raster::extent(x)
    s <- raster::resample(x, s, method='bilinear')
    
    return(s)
}
