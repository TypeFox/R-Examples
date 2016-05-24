## ====================
## Plotting the legend
## ====================

PST.legend <- function(pos, text, colors, cex=1) {

	nbstat <- length(text)

	## Computing some parameters for the legend's plotting
	if (pos=="bottom") {
		if (nbstat > 6) 
			nbcol <- 3
		else 
			nbcol <- 2

		leg.ncol <- ceiling(nbstat/nbcol) 
	}
	else 
		leg.ncol <- 1

	## leg.inset <- -0.2 + ((2-leg.ncol)*0.025)

	## Setting graphical parameters while saving them in savepar
	savepar <- par(mar = c(1, 1, 0.5, 1) + 0.1, xpd=FALSE)

	## Restoring graphical parameters
	on.exit(par(savepar))

	plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
		## legend(position, fill = cpal, legend = ltext, cex = fontsize)

	legend(pos,
		## inset=c(0,leg.inset),
		legend=text,
		fill=colors,
		ncol=leg.ncol,
		bty="o",
		cex=cex)

}
