print.bivden <- function(x, ...){
	cat("Bivariate kernel density/intensity estimate\n\n")
	if(length(unique(x$h))==1) cat("Fixed isotropic smoothing with h =",unique(x$h[!is.na(x$h)]),"unit(s)\n")
	else cat("Adaptive isotropic smoothing with (pilot) h =",x$pilotH,"global h =",x$globalH,"unit(s)\n")
	
	cat("No. of observations:",sum(x$counts),"\n")
}