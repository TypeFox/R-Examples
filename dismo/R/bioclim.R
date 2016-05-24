# Author: Robert J. Hijmans
# Date : December 2009
# Version 1.0
# Licence GPL v3

setClass('Bioclim',
	contains = 'DistModel',
	representation (
		min='vector',
		max='vector'
	),	
	prototype (	
	),
	validity = function(object)	{
		return(TRUE)
	}
)


if (!isGeneric("bioclim")) {
	setGeneric("bioclim", function(x, p, ...)
		standardGeneric("bioclim"))
}	


setMethod('bioclim', signature(x='matrix', p='missing'), 
	function(x, p, ...) {
	
		bioclim(data.frame(x), ...)
		
	}
)

setMethod('bioclim', signature(x='data.frame', p='missing'), 
	function(x, p, ...) {
		for (i in ncol(x):1) {
			if (is.factor(x[,i])) {
				warning('variable "', colnames(x)[i], '" was removed because it is a factor (categorical variable)')
				x <- x[, -i]
			}
		}
		
		bc <- new('Bioclim')

		x = stats::na.omit(x)
		if (ncol(x) == 0) {	stop('no usable variables') 	}
		if (nrow(x) < 2) {	stop('insufficient records') 	}
		
		bc@presence <- x
		bc@min <- apply(x, 2, min)
		bc@max <- apply(x, 2, max)
		bc
	}
)

setMethod('bioclim', signature(x='Raster', p='matrix'), 
	function(x, p, ...) {
		m <- extract(x, p)
		bioclim(data.frame(m), ...)
	}
)

setMethod('bioclim', signature(x='Raster', p='data.frame'), 
	function(x, p, ...) {
		m <- extract(x, p)
		bioclim(data.frame(m), ...)
	}
)

setMethod('bioclim', signature(x='Raster', p='SpatialPoints'), 
	function(x, p, ...) {
		m <- extract(x, coordinates(p))
		bioclim(data.frame(m), ...)
	}
)

setMethod('bioclim', signature(x='SpatialGridDataFrame', p='SpatialPoints'), 
	function(x, p, ...) {
		x <- brick(x)
		p <- coordinates(p)
		bioclim(x, p)
	}
)

setMethod('bioclim', signature(x='SpatialGridDataFrame', p='matrix'), 
	function(x, p, ...) {
		x <- brick(x)
		bioclim(x, p)
	}
)

