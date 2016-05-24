# The 'stacf' function is defined as, per (Pfeifer & Stuart, 1980), the 
# following expression:
# 	acf(l,0,s) = cov(l,0,s) / sqrt( cov(l,l,0) * cov(0,0,0) )

stacf <- function(data, wlist, tlag.max=NULL, plot=TRUE, use.ggplot=TRUE) {

	# If only the weights matrix of first order is specified	
	if (is.matrix(wlist))
		wlist <- list(diag(dim(wlist)[1]), wlist)

	# If no tlag.max is specified
	if (is.null(tlag.max))
		tlag.max <- floor(10 * log10(nrow(data)))
	
	# Call C++ function for optimized computation
	if (is.data.frame(data))
		out <- stacfCPP(as.matrix(data), wlist, tlag.max)
	else
		out <- stacfCPP(data, wlist, tlag.max)

	colnames(out) <- paste("slag", 0:(length(wlist) - 1))
	rownames(out) <- paste("tlag", 1:tlag.max)
	
	# Plot stacf (and still returns the stacf matrix for efficient use)
	if (plot) {
		stplot(out, 2 / sqrt(nrow(data) * ncol(data)), match.call(), ggplot=use.ggplot)
		invisible(out)
	}
	else
		out

}

# A faire :
# - Essayer de separer les fichiers stcov.cpp et stacf.cpp tout en gardant
#   la coherence du stacfCPP (qui utilise la fonction stcovCPP definie dans
#   stcov.cpp).