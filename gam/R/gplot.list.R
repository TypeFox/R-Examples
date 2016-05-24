"gplot.list" <-
function(x, y, se.y = NULL, xlab, ylab, residuals = NULL, rugplot = FALSE, scale = 
	0, se = FALSE, fit = TRUE, ...)
{
	if(length(x) != 2) {
		warning(paste("A perspective plot was requested for \"", ylab,
			"\" but the \"x\" variable has dimension other than 2",
			sep = ""))
		invisible(return(0))
	}
	names(x) <- xlab
	x <- data.matrix(data.frame(x))
	#	UseMethod("gplot")
	gplot.matrix(x, y, se.y, xlab, ylab, residuals, rugplot, scale, se,
		fit, ...)
}
