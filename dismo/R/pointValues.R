# Author: Robert J. Hijmans
# contact: r.hijmans@gmail.com
# Date : December 2009
# Version 0.1
# Licence GPL v3

pointValues <- function(x, p, a, uniquecells=TRUE, na.rm=TRUE) {
	if (inherits(p, 'SpatialPoints')) {
		p <- coordinates(p)
	}
	if (inherits(a, 'SpatialPoints')) {
		a <- coordinates(a)
	}
	pa <- rbind(cbind(1,p), cbind(0,a))
	np <- nrow(pa)
	cell <- cellFromXY(x, pa[,2:3])
	cell <- stats::na.omit(cell)
	if (length(cell) < np) {
		frac <- length(cell) / np
		warning("points on raster=", frac," times total" )
	}
	if (uniquecells) {
		nc <- length(cell)
		cell <- unique(cell)
		if (length(cell) < nc) {
			frac <- length(cell) / nc
			warning("unique cells=", frac," times total" ) 
		}
	} 
	vals <- cbind(pa[,1], extract(x, cell))
	if (is.null(colnames(vals))) {
		colnames(vals) <- c('presabs', 'value')
	} else {
		colnames(vals)[1] <- 'presabs'
	}
	if (na.rm) {
		vals <- stats::na.omit(vals)
		if (length(vals[,1]) < length(cell)) {
			frac <- length(vals[,1]) / length(cell)
			warning("cells with values=", frac," times total unique cells" )
		}
	}
	return(data.frame(vals))
}

