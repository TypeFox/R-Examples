.linIP <- function(ndvi){
	days <- length(ndvi)
	f <- approxfun(x=(-days+1):(2*days), y=c(ndvi,ndvi,ndvi))
	ndvi.interpol <- f(1:days)
	return(ndvi.interpol)
}
