#' Compute a Lorentzian Curve
#' 
#' Computes the y values describing a Lorentzian curve such as seen in an NMR
#' peak. Requires a range of x values and parameters for peak position, area,
#' and gamma (half the peak width at half-height).
#' 
#' 
#' @param x A vector of x values which will be used to compute the
#' corresponding y values.  Use enough to give good resolution.
#'
#' @param x0 The position of the peak.  Must fall in the range of x, of course.
#'
#' @param area The area of the peak, in arbitrary units.
#'
#' @param gamma HWHM, half-width at half-maximum.  The peak "width" in units
#' corresponding to x.
#'
#' @return A vector of y values corresponding to the x values supplied.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @seealso \code{\link{gaussCurve}}, \code{\link{makeSpec}},
#' \code{\link{plotNMRspec}}
#' and \code{\link{plot2DNMRspec}} for drawing NMR spectra.
#'
#' @keywords utilities distributions
#' @export
#' @examples
#' 
#' myx <- seq(0, 100, length.out = 1000) # use lots of point for resolution
#' myy <- lorentzCurve(x = myx, area = 1, x0 = 40, gamma = 5)
#' plot(myx, myy, type = "l", main = "Pure Lorentzian Curve")
#' y = 0.5*max(myy)
#' x = seq(40, 45, 0.5)
#' points(x = x, y = rep(y, length(x)), col = "blue", type = "l")
#' text(x = 42, y = y + 0.005, labels = c("gamma"), col = "blue", srt = 90)
#' 
lorentzCurve <-
function(x, x0, area, gamma) {
	
	# Function to generate Lorentzian curves
	# Bryan Hanson, DePauw Univ, Feb 2011
	# Derived from gaussCurve
	# gamma is peak width @ half height
	# x0 is the center of the peak
	# x is a vector of values @ which the pdf should be computed
	# area is the peak area

	y <- gamma/((x - x0)^2 + gamma^2)
	y <- y*area/pi
	}

