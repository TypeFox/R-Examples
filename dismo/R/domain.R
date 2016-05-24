# Author: Robert J. Hijmans
# contact: r.hijmans@gmail.com
# Date : December 2009
# Version 0.1
# Licence GPL v3



setClass('Domain',
	contains = 'DistModel',
	representation (
		range='vector',
		factors='vector'
	),	
	prototype (	
	),
	validity = function(object)	{
		return(TRUE)
	}
)


if (!isGeneric("domain")) {
	setGeneric("domain", function(x, p, ...)
		standardGeneric("domain"))
}	


setMethod('domain', signature(x='Raster', p='matrix'), 
	function(x, p, ...) {
		m <- extract(x, p)
		domain(as.data.frame(m))
	}
)

setMethod('domain', signature(x='Raster', p='data.frame'), 
	function(x, p, ...) {
		m <- extract(x, p)
		domain(as.data.frame(m))
	}
)

setMethod('domain', signature(x='matrix', p='missing'), 
	function(x, p, ...) {
		domain(as.data.frame(x), ...)
	}
)

setMethod('domain', signature(x='data.frame', p='missing'), 
	function(x, p, factors, ...) {
		
		if (missing(factors)) {
			factors <- colnames(x)[ sapply(x, function(x) is.factor(x)) ]
		}

		x = stats::na.omit(x)
		
		if (ncol(x) == 0) {	stop('no usable variables') 	}
		if (nrow(x) < 2) {	stop('insufficient records') 	}
		
		r <- apply(x, 2, FUN=range)

		d <- new('Domain')
		d@presence <- x
		d@factors <- factors
		d@range <-  abs(r[2,] - r[1,])
		norange = which(d@range == 0)
		if (length(norange) > 0) {
			for (i in length(norange):1) {
				index = norange[i]
				name <- colnames(d@presence)[index]
				if (! name %in% factors) {
					warning('variable "', name, '" was removed because it has no variation for the training points')
					d@presence <- d@presence[, -index]
					d@range <- d@range[-index]
				}
			}
		}
		if (ncol(d@presence) == 0) {
			stop('no usable variables')
		}
		d
	}
)

setMethod('domain', signature(x='Raster', p='SpatialPoints'), 
	function(x, p, ...) {
		m <- extract(x, p)
		domain(as.data.frame(m))
	}
)


