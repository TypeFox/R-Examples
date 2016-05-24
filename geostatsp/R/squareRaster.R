squareRaster = function(x, cells=NULL)  {
	UseMethod("squareRaster")
	
}

squareRaster.matrix = function(x,   cells=NULL) {
	
	x = extent(x)
	squareRaster(x, cells)
	
}
squareRaster.Extent = function(x, cells=NULL) {
	
	if(is.null(cells))
		warning("cells must be specified if x is a matrix or extent")
	cells = as.integer(cells[1])
	x = raster(x, ncol=cells, nrow=cells)
	squareRaster(x, cells)
}

squareRaster.BasicRaster = squareRaster.RasterLayer = function(x, cells=NULL) {
	x=raster(x)
	if(is.null(cells)) {
		cells = ncol(x)
	} else {
		ncol(x) = as.integer(cells[1])
	}
	Ny = ceiling(signif( (ymax(x) - ymin(x))/xres(x), 10) )
	ymax(x) = ymin(x) + Ny*xres(x)
	nrow(x) = Ny
	x
}

squareRaster.SpatialPointsDataFrame = squareRaster.SpatialPoints =
		squareRaster.SpatialPolygons = squareRaster.SpatialPolygonsDataFrame =
		function(x, cells=NULL) {
	
	result = squareRaster(extent(x), cells)
	proj4string(result) = proj4string(x)
	result
}