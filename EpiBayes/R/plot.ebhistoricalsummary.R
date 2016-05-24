#' @title 
#' Plot Method for Summary of EpiBayes Historical Object
#' 
#' @description
#' This function plots summary measurements for posterior distributions of cluster-level 
#'     prevalences across all time periods considered. It does so by plotting the output 
#'     from the \code{\link{summary.ebhistorical}} method, which is of class
#'     \code{ebhistoricalsummary}.
#'
#' @param x An object of class \code{ebhistoricalsummary} (e.g., the output of 
#'     method \code{\link{summary.ebhistorical}}).
#' @param ... Additional arguments to be passed on to plot.
#'
#' @return
#' The summary statistics are plotted in the current plotting window. The colors and 
#'     plotting characters represent the various subzones and the actual points in the plot 
#'     are the summary statistics (vertical axis) measured in a particular time period 
#'     (horizontal axis).
#'
#' @seealso 
#' This is a method for objects of class \code{ebhistoricalsummary} returned by the method
#'     \code{\link{summary.ebhistorical}} and creates its own class of object much like the 
#'     summary method for \code{lm} objects.
#'
#' @export
#' @name plot.ebhistoricalsummary
#' @import scales
plot.ebhistoricalsummary = function(x, ...){
	
		# Derive some other useful values
		time.labels = x$ForOthers$time.labels
		subzones = x$ForOthers$unique.subzones
		n.periods = ncol(x$sumstat.mat)
		n.subzones = nrow(x$sumstat.mat)
		sumstat.mat = x$sumstat.mat
		
	# Plot the summary statistics throughout time
	# Set up vector of default colors
	def.cols = c(
			scales::alpha("black", 0.7), 
			scales::alpha("blue", 0.7), 
			scales::alpha("red", 0.7), 
			scales::alpha("green", 0.7), 
			scales::alpha("magenta", 0.7), 
			scales::alpha("gray50", 0.7), 
			scales::alpha("cyan", 0.7), 
			scales::alpha("orange", 0.7), 
			scales::alpha("forestgreen", 0.7), 
			scales::alpha("blueviolet", 0.7), 
			scales::alpha("gray75", 0.7), 
			scales::alpha("lightblue", 0.7), 
			scales::alpha("indianred1", 0.7), 
			scales::alpha("lightgreen", 0.7), 
			scales::alpha("lightpink", 0.7), 
			rev(topo.colors(n.periods))
			)

	# Declare a new plotting device	
	dev.new()
	
		plot(which(!(sumstat.mat[1, ] %in% c(0,1))), sumstat.mat[1, which(!(sumstat.mat[1, ] %in% c(0,1)))], 
		type = "o", 
		col = scales::alpha("black", 0.7), 
		pch = 15, 
		ylim = c(0, max(sumstat.mat)), 
		xlim = c(1, n.periods),
		main = paste0("Posterior Summary Statistics of Cluster-level Prevalence \n within Subzone over ", n.periods, " Years of Sampling"),
		xlab = "Year",
		ylab = "Cluster-level Prevalence within Subzone",
		xaxt = "n",
		...)
	axis(1, at = 1:n.periods, labels = time.labels)
	
	if(n.subzones > 1){
		for (i in 2:n.subzones){
			lines(which(!(sumstat.mat[i, ] %in% c(0,1))), sumstat.mat[i, which(!(sumstat.mat[i, ] %in% c(0,1)))], type = "o", col = def.cols[i], pch = 15 + i - 1)
		}
	}
	
	legend("topright", legend = subzones, pch = 15:(15+n.subzones), col = def.cols)
	
}