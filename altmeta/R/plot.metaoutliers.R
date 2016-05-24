plot.metaoutliers <- function(x, xtick.cex = 1, ytick.cex = 0.5, ...){
	std.res <- x$std.res
	n <- length(std.res)
	abs.max <- max(abs(std.res))
	abs.limit <- abs.max*1.1
	plot.default(x = std.res, y = n:1, ylab = "Study", xlab = "Standardized Residual",
		ylim = c(1, n), xlim = c(-abs.limit, abs.limit), yaxt = "n", xaxt = "n", pch = 20, cex = 1.5, cex.lab = 1)
	axis(side = 2, at = 1:n, labels = n:1, cex.axis = ytick.cex)
	axis(side = 1, at = -4:4, cex.axis = xtick.cex)
	if(abs.limit > 3) abline(v = c(-3, 3), lwd = 2, lty = 2)
}