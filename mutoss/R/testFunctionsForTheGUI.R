# A File for temporarily functions to test the GUI.
# 
# Author: MarselScheer
###############################################################################


pValuesPlot = function(pValues) {
	plot(ecdf(pValues), do.points=FALSE, verticals = TRUE, main="ecdf", ylim=c(0,1))
	abline(0,1, col=2)
}

adjPValuesPlot = function(adjPValues, alpha) {
	plot(sort(adjPValues), main="Adjusted p-values", ylab="adjusted p-values", xlab="ordered index", ylim=c(0,1))
	if (!missing(alpha)) { 
		abline(alpha,0, col=2) 
	}
}
