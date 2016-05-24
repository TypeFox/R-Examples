plot.klcv <- function(x, ...){
	rho <- x$rho
	out_klcv <- x$klcv
	min_klcv <- x$rhoid
	plot(rho, out_klcv, ...)
	abline(v = rho[min_klcv], lty = 2, lwd = 2, col = 1)
}
