#as.im = function(X, ...) {
#	UseMethod("as.im")
#}

asImRaster = function(X, ...) {
  if(requireNamespace('spatstat',quietly=TRUE)){	
	res=spatstat::as.im(raster::as.matrix(X)[nrow(X):1,], 
			xrange=bbox(X)[1,], 
			yrange=bbox(X)[2,], ...)
} else {
	message("Install the spatstat package if you wish to use this function.")
	res = NULL
}
return(res)
}

#as.im.default = function(X, ...) {
#	spatstat::as.im(X, ...) 
#}