#' @title 
#' Plot Method for EpiBayesHistorical Output
#' 
#' @description
#' This method plots the output of the function \code{\link{EpiBayesHistorical}}. It does 
#'     so by plotting the posterior distributions of the cluster-level prevalence across 
#'     all time periods in different colors within a particular subzone. Similar plots 
#'     are made for each subzone in separate plotting windows (\emph{Caveat: too many 
#'     subzones can make manipulation of the plotting windows unmageable}).
#'
#' @param x An object of class \code{ebhistorical} (e.g., the output of 
#'     function \code{\link{EpiBayesHistorical}}).
#' @param ... Additional arguments to be passed to \code{plot}.
#'
#' @return
#' One plotting window for each subzone including posterior distributions for the 
#'     cluster-level prevalence across all time periods.
#'
#' @export
#' @name plot.ebhistorical
#' @import shape
#' @import scales
plot.ebhistorical = function(x, ...){
	
	# Extract the outputs of the historical object
	periodi.tau = x$RawPost
	periodi.betabust = x$BetaBusterEst
	orig.tauparm = x$ForOthers$orig.tauparm
	subzones = x$ForOthers$unique.subzones
	maxdens = x$ForOthers$maxdens
	burnin = x$ForOthers$burnin
	unique.periods = sort(unique(x$ForOthers$input.df$period))
	n.periods = x$ForOthers$n.periods
	n.subzones = length(x$ForOthers$unique.subzones)

	# Plot the prior and the subsequent posterior distributions with, optionally, the estimated beta distributions 
	# Color indexes the period of the distribution of tau
	# Line type denotes whether the plot is the density estimate (solid) or the BetaBuster fit (dashed) of the distribution of tau
	# Set the palette to topological colors for nice ramp throughout periods
	default.palette = palette()

	for (j in 1:n.subzones){
		
		# Create a new plotting device
		dev.new()
		par(mar = c(5, 4, 4, 6) + 0.1)
		palette(rev(topo.colors(n.periods, 0.7)))
		# Blank plot for the subzone
			# Get the maximum y plotting value and minimum legend plotting value
			maxploty = max(1, maxdens[, j])
			minplotleg = ifelse(maxploty == 1, 0.6, maxploty - 3)
			
		plot(c(0, 0.1), c(0, maxploty), type = "n", main = paste("Inheritance of Information Across \n Periods in Subzone", subzones[j]), xlab = expression(paste("Cluster-level Prevalence (", tau, ")")) , ylab = "Density")
			#Prior
			curve(dbeta(x, orig.tauparm[1], orig.tauparm[2]), from = 0, to = 1, add = TRUE, col = scales::alpha("black", 0.7), lty = 1)
			
			#Posteriors
			for(i in 1:n.periods){
				if (any(periodi.tau[[i]][j, -c(1:burnin)] != 0)){
					lines(density(periodi.tau[[i]][j, -c(1:burnin)], from = 0, to = 1), col = i, lty = 1)
				} 
			}
	
		shape::colorlegend(col = rev(topo.colors(n.periods)), zlim = c(unique.periods[1], unique.periods[n.periods]), zval = unique.periods, posx = c(0.87, 0.90))
		# color.legend(0.08, minplotleg, 0.1, maxploty, legend = c(unique.periods[1], unique.periods[n.periods]), rect.col = rev(topo.colors(n.periods + 1)), gradient = "y")
		# legend("topright", legend = c("Prior", periodi.names), col = c("black", rev(topo.colors(n.periods + 1))), lty = rep(1, n.periods + 1))
	}
	
	# Reset the color palette to the default
	palette(default.palette)

}