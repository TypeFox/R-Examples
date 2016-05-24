# R function for the raster package
# Author: Robert J. Hijmans
# Date : January 2009
# Version 0.9
# Licence GPL v3


cellsFromExtent <- function(object, extent, expand=FALSE) {
	object <- raster(object) 
	extent <- alignExtent(extent(extent), object)
	innerBox <- intersect(extent(object), extent)
	if (is.null(innerBox)) { 
		return(NULL) 
	}
	
	srow <- rowFromY(object, innerBox@ymax - 0.5 * yres(object))
	erow <- rowFromY(object, innerBox@ymin + 0.5 * yres(object))
	scol <- colFromX(object, innerBox@xmin + 0.5 * xres(object))
	ecol <- colFromX(object, innerBox@xmax - 0.5 * xres(object))
	
	if (expand) {
		srow <- srow - round((extent@ymax - innerBox@ymax) / yres(object))
		erow <- erow + round((innerBox@ymin - extent@ymin) / yres(object))
		scol <- scol - round((innerBox@xmin - extent@xmin) / xres(object))
		ecol <- ecol + round((extent@xmax - innerBox@xmax) / xres(object))
	}

	return(cellFromRowColCombine(object, srow:erow, scol:ecol))
}

