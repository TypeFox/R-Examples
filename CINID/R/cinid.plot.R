cinid.plot <- function(cinid.out, breaks = 50, xlab = "Headcapsules width", ylab1 = "Density", ylab2 = "Number of Larvae", main = "", ...) {
  
	if((length(cinid.out) != 2) | (is.list(cinid.out) == FALSE))
		stop("Erreur : cinid.out must be a list of two elements provided by the cinid.table function")

	Mindiv <- cinid.out$indiv
	Mpop <- cinid.out$pop
  par(mar = c(5, 6, 3, 5))
	H <- hist(as.numeric(as.vector(Mindiv$HCW)), breaks = 50, plot = F)

	## define the four unimodal distribution
	X <- seq(min(H$breaks), max(H$breaks), length.out = 1000)
	mod1 <- (exp(- 0.5 * ((X - Mpop[1, 1]) / Mpop[2, 1]) ^ 2) / (Mpop[2, 1] * sqrt(2 * pi)))
	mod2 <- (exp(- 0.5 * ((X - Mpop[1, 2]) / Mpop[2, 2]) ^ 2) / (Mpop[2, 2] * sqrt(2 * pi)))
	mod3 <- (exp(- 0.5 * ((X - Mpop[1, 3]) / Mpop[2, 3]) ^ 2) / (Mpop[2, 3] * sqrt(2 * pi)))
	mod4 <- (exp(- 0.5 * ((X - Mpop[1, 4]) / Mpop[2, 4]) ^ 2) / (Mpop[2, 4] * sqrt(2 * pi)))

	## maximal density
	maxi <- max(max(H$density), max(Mpop[6, 1] * mod1), max(Mpop[6, 2] * mod2), max(Mpop[6, 3] * mod3), max(Mpop[6, 4] * mod4))

	## plot histogram
	H <- hist(as.numeric(as.vector(Mindiv$HCW)), breaks = breaks, plot = T, prob = T, xlab = xlab, ylab = "", main = main, xaxt = "n", yaxt = "n", ylim = c(0, maxi))
	axis(1, at = seq(min(H$breaks), max(H$breaks), by = (H$breaks[2] - H$breaks[1]) * 2), cex.axis = 0.7, line = -0.5)
	yaxis <- axis(2, las = 1)
	mtext(side = 2, text = ylab1, line = 4)
	axis(4, yaxis, labels = pretty(0:max(H$counts), 5)[1:length(yaxis)], las = 1)
	mtext(side = 4, text = ylab2, line = 3)

	## color individual with intermined instar
	par(new = T)
	indet <- Mindiv[which(Mindiv$instar_determ == "Indet"), 1]
	h <- hist(as.numeric(as.vector(indet)), breaks = H$breaks, plot = F)
	v <- h$counts / sum(H$counts) / (H$breaks[2] - H$breaks[1])
	rect(H$mids - 0.5 * (H$breaks[2] - H$breaks[1]), 0, H$mids + 0.5 * (H$breaks[2] - H$breaks[1]), v, col = "grey")

	## plot four probabilities' densities
	lines(X, Mpop[6, 1] * mod1, lwd = 1.5)
	lines(X, Mpop[6, 2] * mod2, lwd = 1.5)
	lines(X, Mpop[6, 3] * mod3, lwd = 1.5)
	lines(X, Mpop[6, 4] * mod4, lwd = 1.5)
	invisible(H)
}
