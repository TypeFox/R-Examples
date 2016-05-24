
# ratesHistogram <- function(...)
		# phylorates = a plot.bammdata object
		# plotBrks = boolean, should breaks be plotted
		# xlab = x-axis label
		# ylab = y-axis label
		# lwd = lwd for breaks lines
		# lty = lty for breaks lines
		# brksCol = color for breaks lines
		# ... additional arguments passed on to mtext for axis labels

ratesHistogram <- function(phylorates, plotBrks = TRUE, xlab = 'speciation rate', ylab = 'density', lwd = 0.2, lty = 1, brksCol = 'black', ...) {
	
	if (!identical(names(phylorates), c("coords", "colorbreaks", "palette", "colordens"))) {
		stop("phylorates must be a saved plot.bammdata object.")
	}

	plot.new();
	x <- phylorates$colordens[,1];
	y <- phylorates$colordens[,2];
	plot.window(xlim = c(min(x), max(x)), ylim = c(0, max(y)));
	segments(x, y, x, 0, lend = 2, col = phylorates$colordens[,3], lwd = 3);
	axis(1, signif(seq(min(0, min(x)), max(x), length.out = 5), 2), xaxs = "i", cex.axis = 0.75, tcl = NA, mgp = c(0, 0.25, 0));
	axis(2, round(c(-1, seq(0, max(y), length.out = 4))), las = 1, yaxs = "i", cex.axis = 0.75, tcl = NA, mgp = c(0, 0.35, 0));
	
	mtext(xlab, side = 1, line = 1.5, ...);
	mtext(ylab, side = 2, line = 1.5, ...);
	
	#add breaks as vertical lines
	if (plotBrks) {
		abline(v = phylorates$colorbreaks, lwd = lwd, lty = lty, col = brksCol);
	}
}





