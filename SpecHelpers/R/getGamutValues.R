##' Look up gamut and white point values in the 1931 CIE system
##'
##' These functions provide a simple way of storing white point and
##' gamut data for use in drawing CIE chromaticity diagrams.
##' 
##' @param gamut A character string giving the name of the desired gamut.  One of
##' \code{c("Apple", "CIE", "Adobe", "sRGB", "NTSC", "SWOP")}.
##'
##' @param white The desired white point value.  One of \code{c("D65", "E", "C", "D50")}.
##'
##' @seealso \code{\link{plotCIEchrom}} for examples of this function in use.
##' 
##' @return A data frame with columns x, y containing the vertices of the
##' requested gamut in CIE chromaticity coordinates, or, for a white point,
##' a data frame containing the coordinates of the requested white point.
##'
##' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
##'
##'
##' @aliases getGamutValues getWhiteValues
##'
##' @export
##' @keywords utilities
##'


getGamutValues <- function(gamut) {
	
	# gamut data verified with multiple sources
	# CIE is the theoretical primaries from CIE
	# sRGB is the standard for computer monitors and HDTV
	# Apple is the standard for Apple monitors based upon Trinitron phosphors
	#       which is no longer in use. Same as SGI.
	# NTSC is the 1953 CRT TV standard

	# Bryan Hanson, DePauw University, March 2013 hanson@depauw.edu
	
	# Most or all of these gamuts corners are in order of r, g, b (ccw)
	
	if (gamut == "Apple") { # relative to D65
		x <- c(0.640, 0.300, 0.150)
		y <- c(0.330, 0.600, 0.060)
		ans <- data.frame(x, y)
		}
		
	if (gamut == "CIE") { # relative to E
		x <- c(0.735, 0.274, 0.167)
		y <- c(0.265, 0.717, 0.009)
		ans <- data.frame(x, y)
		}

	if (gamut == "Adobe") { # relative to D65
		x <- c(0.640, 0.210, 0.150)
		y <- c(0.330, 0.600, 0.060)
		ans <- data.frame(x, y)
		}

	if (gamut == "sRGB") { # relative to D65
		x <- c(0.640, 0.300, 0.150)
		y <- c(0.330, 0.600, 0.060)
		ans <- data.frame(x, y)
		}

	if (gamut == "NTSC") { # relative to C
		x <- c(0.670, 0.210, 0.140)
		y <- c(0.330, 0.710, 0.080)
		ans <- data.frame(x, y)
		}

	if (gamut == "SWOP") { # A type of CMYK (approx values)
		x <- c(0.205, 0.172, 0.225, 0.430, 0.610, 0.470)
		y <- c(0.125, 0.226, 0.540, 0.500, 0.320, 0.235)
		ans <- data.frame(x, y)
		}

	ans
	}
	
##' @rdname getGamutValues

getWhiteValues <- function(white) {
	
	# Bryan Hanson, DePauw University, March 2013 hanson@depauw.edu

	if (white == "D65") ans <- data.frame(x = 0.3127, y = 0.3290)
	if (white == "E") ans <- data.frame(x = 0.333, y = 0.333)
	if (white == "C") ans <- data.frame(x = 0.3101, y = 0.3161)
	if (white == "D50") ans <- data.frame(x = 0.3457, y = 0.3585)

	ans
	}