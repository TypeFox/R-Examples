plot.sglasso <- function(x, ...){
	rho <- x$rho
	nv <- x$nv
	theta <- x$theta[-(1:nv), , drop = FALSE]
	grd <- x$grd[-(1:nv), , drop = FALSE]
	w <- x$w
	grd <- grd / w
	matplot(rho, t(theta), col = 1, type = "l", lty = 1, xlab = expression(rho), ylab = "", main = "Coefficients Path")
	op <- par(ask = dev.interactive())
	matplot(rho, abs(t(grd)), col = 1, type = "l", lty = 1, xlab = expression(rho), ylab = "", main = "Weighted Score Path")
	par(op)
}
