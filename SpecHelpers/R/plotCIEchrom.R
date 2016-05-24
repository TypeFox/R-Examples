
##' Draw the 1931 CIE chromaticity diagram
##'
##' This function draws the 1931 CIE chromaticity diagram with
##' various decorations and annotations.  
##' 
##' @param gradient Character: either \code{"sl"}, \code{NULL}, or a
##' data frame with columns x and y.  If \code{NULL}, no gradient is drawn.
##' If \code{"sl"} a gradient filling the entire spectral locus is drawn.
##' If a data frame, the vertices should specify a polygon to be filled with
##' the gradient (see the examples for convenient ways to specify the gradient).
##' 
##' @param colSpace Character string giving the color space to use for
##' drawing the gradient.  One of \code{c("sRGB", "Apple RGB")}.
##' \code{Apple RGB} is mainly of historical interest; no physical
##' devices use it at this time.
##' 
##' @param  ex Numeric.  The 'exposure' to use.  The exposure must be
##' used with \strong{extreme care}.  Larger values of \code{exposure}
##' make the white point whiter in the plot, and lightens colors near
##' the spectral locus (driving some off the plot!).  The purpose is to
##' alter the aesthetics of the plot - that is, to make the white "whiter"
##' so that it looks "right".  The effect of exposure will vary with the display device.
##' 
##' @param opts A character vector of options to be employed.  One or
##' more of c("D65", "D50", "C", "E", "specLocus", "purples", "Munsell",
##' "sRGB", "SWOP", "Apple", "NTSC", "Adobe", "CIE").  The first few of
##' these are reference white points.  \code{"specLocus"} and \code{"purples"}
##' cause the spectral locus and line of purples to be labeled.  \code{"Munsell"}
##' causes the approximate Munsell hues to be marked along the spectral
##' locus at the appropriate wavelength.  The last few options cause the requested gamut to be outlined.
##' 
##' @param title A character string to be plotted at the top of the diagram.
##' If NULL, the title defaults to "1931 CIE Chromaticity Diagram".
##' If no title is desired, set it to an empty string.
##' 
##' @param \dots Additional arguments to be passed downstream, to \code{grid} functions.
##' 
##' @return A plot is drawn using \code{grid} graphics.
##' 
##' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
##'
##' @references For \code{opts = "Munsell"}
##' the Munsell designation by wavelength are taken from Romney & Indow
##' \url{pnas.org/cgi/doi/10.1073/pnas.162368999}
##'
##' @section Warning:
##' The appearance of the color gradient will vary with the
##' device, surface and incident light used to view it and is not likely correct
##' anywhere.  \strong{The appearance varies strongly with exposure}.
##'
##' @export
##' @keywords hplot
##'
##' @importFrom grid arrow grid.points grid.polygon unit grid.raster grid.xaxis grid.yaxis
##' @importFrom grid pushViewport viewport grid.text gpar grid.newpage
##' @importFrom grDevices as.raster
##' @importFrom utils data
##'
##' @examples
##' require("grid")
##' plotCIEchrom() # no gradient
##' ## These are a too slow for CRAN checks:
##' \dontrun{
##' plotCIEchrom(gradient = "sl") # basic plot
##' # Notice there is not much yellow in that plot.  Increase
##' # the exposure to bring in some yellow, at the expense of some blues:
##' plotCIEchrom(gradient = "sl", ex = 1.4)
##' # Next show a gradient for the CMYK printing process
##' # and outline the colors a typical monitor can display.
##' plotCIEchrom(gradient = getGamutValues("SWOP"), opts = c("D65", "SWOP", "sRGB"))
##' }

plotCIEchrom <- function(gradient = NULL, colSpace = "sRGB", ex = 1.0,
	opts = c("D65", "specLocus", "purples"), title = NULL, ...) {

	# Function to draw the CIE chromaticity diagram
	# with various decorations
	
	# Bryan Hanson, DePauw University, March 2013 hanson@depauw.edu
	
##### Get and prepare the data
	# These are the coordinates of the spectral locus, which is a curve
	# describing the pure colors of the spectrum/rainbow

	# Note z = 1- x - y
	# Cutoff the data at 650 nm; beyond that the curve strangely
	# turns back on itself
	keep <- which(CIExyz$wavelength <= 650) # Note CIExyz is lazy loaded
	Lxyz <- CIExyz[keep,]
		
##### Prepare the raster with the color gradient
	
	if (!is.null(gradient)) {
		if (!(colSpace == "sRGB") || (colSpace == "Apple RGB")) stop("colSpace must be sRGB or Apple RGB")
		if (is.character(gradient)) { # captures, crudely, gradient = "sl"
			finras <- as.raster(prepCIEgradient(vertices = Lxyz,
				colSpace = colSpace, ex = ex, ...))
			}
			
		if (is.data.frame(gradient)) {
			finras <- as.raster(prepCIEgradient(vertices = gradient,
				colSpace = colSpace, ex = ex, ...))
			}
		}
	
##### Create the plot using grid graphics

	off <- 0.1 # needed due to vp origin at -0.1
	Lxyz$x <- Lxyz$x + off
	Lxyz$y <- Lxyz$y + off

	# First plot titles & labels in the vp of the entire device

	if (is.null(title)) title <- "1931 CIE Chromaticity Diagram"
	grid.newpage()
	grid.text(title, x = 0.5, y = 0.9,
		gp = gpar(fontface = "bold", cex = 1.2))
	grid.text(expression(italic(x)), x = 0.5, y = 0.05)
	grid.text(expression(italic(y)), x = 0.05, y = 0.5, rot = 90)

	# Now the data in it's own viewport; raster goes underneath (first)

	pushViewport(viewport(width = 0.7, height = 0.7,
		xscale = c(-0.1, 0.9), yscale = c(-0.1, 0.9)))

	if (!is.null(gradient)) {
		grid.raster(finras, x = 0.5, y = 0.5, interpolate = FALSE, default.units = "npc")
		msg1 <- "Warning: the color gradient appearance\nwill vary with the device, surface\n& incident light used to view it\nand is not likely correct anywhere"
		grid.text(msg1, x = 0.98, y = 0.9, gp = gpar(fontface = "italic", cex = 0.9),
			just = "right")
		msg2 <- paste("ex =", ex, "    color space =", colSpace, sep = " ")
		grid.text(msg2, x = 0.98, y = 0.02, gp = gpar(cex = 0.5),
			just = "right")
		}

	grid.polygon(Lxyz$x, Lxyz$y)
	grid.rect()
	tickpos <- seq(0.0, 0.8, by = 0.1)
	grid.xaxis(at = tickpos)
	grid.yaxis(at = tickpos)

	# Labels for the spectral locus
	
	specl <- c(11, 86, 111, 136, 161, 186, 211, 236, 261)
	labs <- c("400 nm  ", "475 nm  ", "500 nm  ", "  525 nm", "  550 nm",
		"  575 nm", "  600 nm", "  625 nm", "  650 nm")
	
	grid.points(x = Lxyz$x[specl], y = Lxyz$y[specl], gp = gpar(col = "black"),
		size = unit(0.5, "char"), default.units = "npc")
	grid.text(label = labs, Lxyz$x[specl], Lxyz$y[specl],
		hjust = c(1, 1, 1, 0, 0, 0, 0, 0, 0),
		vjust = c(1, 0, 0, 0, 0, 0, 0, 0, 1),
		gp = gpar(cex = 0.75))
	
	# Optional labeling
	
	# Need to add lty and a legend to the various gamut options
	
	# White point labels
	
	if ("D65" %in% opts) {
		wh <- getWhiteValues("D65")
		grid.points(wh$x, wh$y, gp = gpar(col = "black"),
			size = unit(0.5, "char"), default.units = "native")
		grid.text(wh$x, wh$y, label = "  D65", just = "left",
			gp = gpar(cex = 0.75), default.units = "native")
		}

	if ("D50" %in% opts) {
		wh <- getWhiteValues("D50")
		grid.points(wh$x, wh$y, gp = gpar(col = "black"),
			size = unit(0.5, "char"), default.units = "native")
		grid.text(wh$x, wh$y, label = "  D50", just = "left",
			gp = gpar(cex = 0.75), default.units = "native")
		}

	if ("C" %in% opts) {
		wh <- getWhiteValues("C")
		grid.points(wh$x, wh$y, gp = gpar(col = "black"),
			size = unit(0.5, "char"), default.units = "native")
		grid.text(wh$x, wh$y, label = "  C", just = "left",
			gp = gpar(cex = 0.75), default.units = "native")
		}

	if ("E" %in% opts) {
		wh <- getWhiteValues("E")
		grid.points(wh$x, wh$y, gp = gpar(col = "black"),
			size = unit(0.5, "char"), default.units = "native")
		grid.text(wh$x, wh$y, label = "  E", just = "left",
			gp = gpar(cex = 0.75), default.units = "native")
		}

	# gamut outlining
	
	if ("sRGB" %in% opts) {
		g <- getGamutValues("sRGB")
		grid.polygon(g$x, g$y, default.units = "native")
		}

	if ("SWOP" %in% opts) {
		g <- getGamutValues("SWOP")
		grid.polygon(g$x, g$y, default.units = "native")
		}

	if ("Apple" %in% opts) {
		g <- getGamutValues("Apple")
		grid.polygon(g$x, g$y, default.units = "native")
		}

	if ("NTSC" %in% opts) {
		g <- getGamutValues("NTSC")
		grid.polygon(g$x, g$y, default.units = "native")
		}

	if ("Adobe" %in% opts) {
		g <- getGamutValues("Adobe")
		grid.polygon(g$x, g$y, default.units = "native")
		}

	if ("CIE" %in% opts) {
		g <- getGamutValues("CIE")
		grid.polygon(g$x, g$y, default.units = "native")
		}

	# Misc decorations
	
	if ("Munsell" %in% opts) {
		mun <- c(81, 99, 112, 123, 166, 191, 211, 266)
		labs <- c("PB  ", "B  ", "BG  ", "G", "  GY",
			"  Y", "YR", "R")
		
		# grid.points(x = Lxyz$x[sl], y = Lxyz$y[sl], gp = gpar(col = "black"),
			# size = unit(0.5, "char"), default.units = "npc")
		grid.text(label = labs, Lxyz$x[mun], Lxyz$y[mun],
			hjust = c(1, 1, 1, 1, 0, 0, 0.5, 0),
			vjust = c(1, 1, 0, -1, 0, 0, -1.5, 2),
			gp = gpar(cex = 0.75))
		}

	if ("specLocus" %in% opts) {
		grid.text(0.75, 0.55, label = "spectral\nlocus", just = "left",
			gp = gpar(cex = 0.75), default.units = "native")
		grid.segments(0.73, 0.53, 0.61, 0.41, default.units = "native",
			arrow = arrow(ends = "last", length = unit(0.025, "npc"),
			angle = 15, type = "closed"))
		}


	if ("purples" %in% opts) {
		grid.text(0.65, 0.05, label = "line of\npurples", just = "left",
			gp = gpar(cex = 0.75), default.units = "native")
		grid.segments(0.63, 0.07, 0.5, 0.17, default.units = "native",
			arrow = arrow(ends = "last", length = unit(0.025, "npc"),
			angle = 15, type = "closed"))
		}

	if (is.null(gradient)) ans <- NULL
	if (!is.null(gradient)) ans <- finras
	invisible(ans)
	}