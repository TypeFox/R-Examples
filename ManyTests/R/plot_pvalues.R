plot_pvalues <-
function(p) {
	n <- length(p)
	o <- ordered_values(n)
	plot(o, -log(p[order(p, decreasing = TRUE)]), xlab = "Expected values", ylab = "-log(p)")
	}
