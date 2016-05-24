as.bathy <- function(x){

	if (is(x,"bathy")) stop("Object is already of class 'bathy'")

	if (is(x,"SpatialGridDataFrame")) x <- raster::raster(x)

	# if x is a RasterLayer do this
	if (is(x,"RasterLayer")) {
		lat.min <- x@extent@xmin
		lat.max <- x@extent@xmax
		lon.min <- x@extent@ymin
		lon.max <- x@extent@ymax
		
		nlat <- x@ncols
		nlon <- x@nrows
		
		lon <- seq(lon.min, lon.max, length.out = nlon)
		lat <- seq(lat.min, lat.max, length.out = nlat)
		
		bathy <- t(raster::as.matrix(raster::flip(x,direction="y")))
		colnames(bathy) <- lon
		rownames(bathy) <- lat
	}
	
	# if not, it has to be a 3-column table (xyz format)
	if (ncol(x)==3 & !exists("bathy")) {
		bath <- x
	    bath <- bath[order(bath[, 2], bath[, 1], decreasing = FALSE), ]

	    lat <- unique(bath[, 2]) ; bcol <- length(lat)
	    lon <- unique(bath[, 1]) ; brow <- length(lon)

		if ((bcol*brow) == nrow(bath)) {
			bathy <- matrix(bath[, 3], nrow = brow, ncol = bcol, byrow = FALSE, dimnames = list(lon, lat))
			} else {
				colnames(bath) <- paste("V",1:3,sep="")
				bathy <- reshape2::acast(bath, V1~V2, value.var="V3")
			}
	}

	if (!exists("bathy")) stop("as.bathy requires a 3-column table, or an object of class RasterLayer or SpatialDataFrame")

	ordered.mat <- check.bathy(bathy)
	class(ordered.mat) <- "bathy"
	return(ordered.mat)

}
