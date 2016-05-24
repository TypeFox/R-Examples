`cusp.nc.vec` <- 
function (alpha, beta, ..., keep.order = FALSE) 
{
	if(is.loaded('cuspnc')){
		cusp.nc.c(alpha, beta, ..., keep.order=keep.order)
	}
	else{
		cusp.nc(alpha, beta, ...)
	}
}