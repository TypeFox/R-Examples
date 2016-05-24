# The 'stpacf' function is computed by, per (Pfeifer & Stuart, 1980), the 
# Durbin method, i.e. iteratively computed the solution to Yule-Walker's
# system for spatial and time orders growing, and keeping the last parameter.

stpacf <- function(data, wlist, tlag.max=NULL, plot=TRUE, use.ggplot=TRUE) {

	if (is.matrix(wlist))
		wlist <- list(diag(dim(wlist)[1]), wlist)

	# If no tlag.max is specified
	if (is.null(tlag.max))
		tlag.max <- floor(10 * log10(nrow(data)))
	
	if (is.data.frame(data))
		out <- stpacfCPP(as.matrix(data), wlist, tlag.max)
	else
		out <- stpacfCPP(data, wlist, tlag.max)

	colnames(out) <- paste("slag", 0:(length(wlist) - 1))
	rownames(out) <- paste("tlag", 1:tlag.max)

	# Plot stpacf (and still returns the stacf matrix for efficient use)
	if (plot) {
		stplot(out, 2 / sqrt(nrow(data) * ncol(data)), match.call(), ggplot=use.ggplot)
		invisible(out)
	}
	else
		out

}

# To do:
# - Pas oblige d'exporter les fonctions intermediaires stmatCPP_, stmatCPP
#   et stvecCPP.