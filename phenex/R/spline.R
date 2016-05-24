.spline <- function(ndvi){
	days <- length(ndvi)

	f <- splinefun(x=(-days+1):(2*days), y=c(ndvi,ndvi,ndvi), method="monoH.FC")
	ndvi.interpol <- f(1:days)

	ndvi.interpol[which(ndvi.interpol < 0)] <- 0
	ndvi.interpol[which(ndvi.interpol > 1)] <- 1

	return(ndvi.interpol)
}
