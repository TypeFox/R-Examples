# Author: Robert J. Hijmans
# contact: r.hijmans@gmail.com
# Date : December 2009
# Version 0.1
# Licence GPL v3


if (!isGeneric("plot")) {
	setGeneric("plot", function(x, y, ...)
		standardGeneric("plot"))
}	


setMethod("plot", signature(x='DistModel', y='numeric'), 
	function(x, y=1, ...) {
		x <- x@presence
		plot(sort(x[,y]), ylab=colnames(x@presence[y]), ...)
	}
)

if (!isGeneric("points")) {
	setGeneric("points", function(x, ...)
		standardGeneric("points"))
}	

setMethod("points", signature(x='DistModel'), 
	function(x, y=2, ...) {
		x <- x@presence
		points(sort(x[,y]), ...)
	}
)
