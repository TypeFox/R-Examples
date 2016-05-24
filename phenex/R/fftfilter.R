.fftfilter <-
function(ndvi, filter.threshold=3){
	if (!is.numeric(filter.threshold)){ stop("'filter.threshold' should be of type 'numeric'") }
	days <- length(ndvi)
	ndvi <- .linIP(ndvi)

	#Fast Fourier Transfusion
	ndvi.fft <- fft(ndvi)

	#Filter
	ndvi.filtered <- ifelse(abs(Re(ndvi.fft)) < filter.threshold, 0, ndvi.fft)

	#inverse
	ndvi.inv <- fft(ndvi.filtered, inverse=TRUE)
	ndvi.inv.real <- vector(mode="numeric", length=days)
	for (i in 1:days){
		ndvi.inv.real[i] <- (1/days)*sqrt(Re(ndvi.inv[i])^2+Im(ndvi.inv[i])^2)
	}

	model <- ndvi.inv.real

	return(model)
}
