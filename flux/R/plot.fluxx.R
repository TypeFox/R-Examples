plot.fluxx <- 
function(x, ...){
	nmes <- names(x)
	x <- x[[1]]
	plot(ghg ~ time, data = x$inn, ...)
	mtext(nmes, line=0, cex=0.9)
	points(ghg  ~ time, data = x$inn[x$inn$hff,], pch=20, col="red4", cex=0.6)
	points(ghg  ~ time, data = x$mod$model, pch=20, col="blue4", cex=0.7)
	abline(x$mod, col="red4")
	lp <- ifelse(coef(x$mod)[2] < 0, "topright", "topleft")
	with(x$fluss, legend(lp, legend=c(paste("flux =", round(flux,3), unit, "/ m2*h"), paste("R2 =", round(r2,3)), paste("podpu =", round(podpu,2))), bty="n", cex=1.1))
}