replot <- function(part2.list, xlim = NULL, ylim = c(0, 1),
	main = "Curve", legend.location = "topleft", bucket.lim = c(0, .9)) {
	buck            <- part2.list$significant
	time.all        <- part2.list$time.all
	curve.g1 <- part2.list$curve.g1
	curve.g2 <- part2.list$curve.g2
	curve.sd1       <- part2.list$curve.sd1
	curve.sd2       <- part2.list$curve.sd2
	alpha           <- part2.list$alpha["alpha"]
	N.g1            <- part2.list$N.g1
	N.g2            <- part2.list$N.g2
	sig             <- part2.list$sig
	groups          <- part2.list$groups
	
	if(is.null(xlim)) xlim <- c(min(part2.list$time.all), max(part2.list$time.all))
	
	ticks <- seq(xlim[1], xlim[2], round(diff(xlim) / 10 + xlim[1]))
	plot(NULL, ,xlim = xlim, ylim = ylim, ylab = 'Proportion of Fixations',
		xlab = 'Time', axes = FALSE, main = main)
	axis(1, at = ticks)
	axis(2)
	box()
	legend(legend.location, lty = 1:2, legend = groups)
	
	#Make significant area yellow
	buck <- bucket(sig, time.all, ylim = bucket.lim)
	
	#Plot overall estimate of curves
	lines(time.all, curve.g1, lty = 1, lwd = 2)
	lines(time.all, curve.g2, lty = 2, lwd = 2)
	
	#Plot Confidence Interval for Group 1 curve
	lines(time.all, curve.g1 - curve.sd1 * qt(alpha / 2, N.g1 - 1), lty = 1, lwd = 1,
		col = "gray44")
	lines(time.all, curve.g1 + curve.sd1 * qt(alpha / 2, N.g1 - 1), lty = 1, lwd = 1,
		col = "gray44")
	
	#Plot Confidence Interval for Group 2 curve
	lines(time.all, curve.g2 - curve.sd2 * qt(alpha / 2, N.g2 - 1), lty = 2, lwd = 1,
		col = "gray44")
	lines(time.all, curve.g2 + curve.sd2 * qt(alpha / 2, N.g2 - 1), lty = 2, lwd = 1,
		col = "gray44")
}