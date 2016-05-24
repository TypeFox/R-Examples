
#' @title
#' Plot a pretty electromagnetic spectrum
#'
#' @description
#' This function plots an annotated electromagnetic spectrum.
#' There are options to include annotations about the molecular
#' effects and/or typical applications in technology.
#'
#' @param molecular Logical. Add annotations about molecular effects?
#' 
#' @param applications Logical. Add annotations about applications?
#' 
#' 
#' @return None. Side effect is a plot.
#'
#' @section Note:
#' The diagram is wider than a standard \code{R} graphics device.
#' You should send it to a \code{pdf} or similar device
#' with the width set to 11" or so.
#' 
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @export
#'
#' @importFrom grid grid.newpage viewport upViewport pushViewport
#' grid.rect gpar grid.lines grid.xaxis grid.text grid.segments
#'
#' @section Details:
#' Obviously not to scale, but hopefully aesthetically pleasing!
#' 
#' @keywords utilities
#' 
#' @examples
#' dev.new(width = 10.5, height = 3)
#' emSpectrum()
#' emSpectrum(molecular = FALSE, applications = FALSE)
#' dev.off()

emSpectrum <- function(molecular = TRUE, applications = TRUE) {

	if (!requireNamespace("grid", quietly = TRUE)) {
		stop("You need to install package grid to use this function")
		}	

	n <- 8 # no. of cycles to be strung together
	x <- seq(0, n*2*pi, length.out = n*360) # one cycle described by 360 data points
	x2 <- x*exp(x/20) # Stretch the x-axis
	
	# Establish the main viewport
	# Viewport sizes in inches, but positions in npc
	grid.newpage()
	mainVP <- viewport(width = unit(10, "inches"), height = unit(3, "inches"))
	pushViewport(mainVP)


	# Draw the visible spectrum as a series of rectangles
	# in the main viewport
	# These need to be underneath the other items, so draw first
	le <- 0.52 # left edge of spectrum
	w <- 0.007 # width of a particular rectangle
	speccols <- c("red", "orange", "yellow", "green", "blue", "darkviolet")
	for (i in 1:length(speccols)) {
		grid.rect(x = (le + w*(i-1)), y = 0.65,
		width = w, height = 0.2,
		gp = gpar(fill = speccols[i], col = "transparent"))
		}
	
	# Overlay the data in its own viewport
	specVP <- viewport(width = unit(10, "inches"), height = unit(0.5, "inches"),
		y = 0.65)
	pushViewport(specVP)
	
	xx <- 1 - x2/max(x2) # moves long wavelengths to the left side
	yy <- 0.5*sin(x) + 0.5
	grid.lines(xx, yy, gp = gpar(lwd = 1.5))
	
	# Draw the axis in its own viewport
	# Labeling the ticks based on http://stackoverflow.com/a/33133919/633251
	# Thanks Josh.
	upViewport() # back to mainVP
	axisVP <- viewport(width = unit(10, "inches"), height = unit(0.1, "inches"),
		y = 0.65)
	pushViewport(axisVP)
	## First plot an axis with all ticks and no labels
	ticks_at <- seq(1/15, 14/15, length.out = 13)
	grid.xaxis(at = ticks_at, label = FALSE)
	grid.lines(x = c(0, 1/5), y = c(0, 0))
	grid.lines(x = c(14/15, 1), y = c(0, 0))
	
	## Then add tick labels where you want them
	labs <- parse(text=c("1", "10^{-3}", "10^{-6}", "10^{-9}", "10^{-12}"))
	labs_at <- seq(1/15, 14/15, length.out = 5)
	grid.xaxis(at = labs_at, label = labs, gp = gpar(cex = 1.0))
	
	# Annotations in the main viewport
	upViewport() # back to mainVP
	grid.text(label = "The Electromagnetic Spectrum", 
		x = 0.5, y = 0.95,
		gp = gpar(cex = 1.2, fontface = "bold"))
	grid.text(label = c("meters", "millimeters", "micrometers", "nanometers", "picometers"),
		x = seq(1/15, 14/15, length.out = 5), y = 0.48)
	grid.text(label = c("radio waves", "microwaves", "infrared", "visible", "ultraviolet", "x-rays", "gamma rays"), gp = gpar(fontface = "italic"),
		x = c(0.08, 0.18, 0.42, 0.54, 0.61, 0.78, 0.93), y = 0.8)
	grid.segments(x0 = 0.03, x1 = 0.0, y0 = 0.8, y1 = 0.8,
		arrow = arrow(angle = 10, length = unit(0.1, "inches")), gp = gpar(lwd = 1.5))
	grid.segments(x0 = 0.9, x1 = 1.0, y0 = 0.9, y1 = 0.9, default.units = "npc",
		arrow = arrow(angle = 5, length = unit(0.25, "inches")), gp = gpar(lwd = 1.5))
	grid.text(label = expression({E==h*nu}==over(hc, lambda)), x = 0.83, y = 0.9, gp = gpar(cex = 1.2))

	# If both molecular and applications are TRUE, we need 2 vp's for annotations.
	# If only one is true, we need just one vp.

	two <- ((molecular) & (applications))
	
	upperVP <- viewport(width = unit(10, "inches"), height = unit(0.5, "inches"),
	y = 0.35)
	lowerVP <- viewport(width = unit(10, "inches"), height = unit(0.5, "inches"),
	y = 0.15)	
	
	if (molecular) {
		pushViewport(upperVP)
		grid.rect(, width = 0.98, gp = gpar(col = "grey"))
		grid.text(label = c("spin transitions", "tumbling", "bond & angle vibration",
			expression(paste(pi %->% pi,"*",~transitions)), "random bond breaking"),
			x = c(0.077, 0.2, 0.42, 0.6, 0.9), y = 0.33)
		grid.text(label = c("nuclear", "electron spin transitions"), x = c(0.05, 0.2), y = 0.66)
		}
	
	if (applications) {
		if (two) {upViewport(); pushViewport(lowerVP)}
		if (!two) pushViewport(upperVP)
		grid.rect(, width = 0.98, gp = gpar(col = "gray"))
		grid.text(label = c("cell", "TV", "wi-fi", "ESR", "motion detector", "insect vision",
			"medical imaging", "decay"),
			x = c(0.04, 0.08, 0.12, 0.2, 0.4, 0.64, 0.8, 0.94), y = 0.33)
		grid.text(label = c("FM", "NMR", "MRI", "microwave oven", "night vision", "human vision",
			"skin cancer", "crystallography", "radioactive"),
			x = c(0.03, 0.07, 0.11, 0.2, 0.4, 0.54, 0.64, 0.8, 0.94), y = 0.66)
		}
	} # end of emSpectrum
	

	