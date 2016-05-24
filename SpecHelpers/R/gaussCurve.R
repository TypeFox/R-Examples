#' Compute a Gaussian Curve
#' 
#' Computes the y values describing a Gaussian distribution given a range of x
#' values and parameters for mu, sigma, and area.  A tail may be introduced
#' into the curve to simulate the behavior of some chromatography peaks.
#' 
#' @param x A vector of x values which will be used to compute the
#' corresponding y values.  Use enough to give good resolution.
#'
#' @param area The area of the peak, in arbitrary units.
#'
#' @param mu The position of the peak.  Must fall in the range of x, of course.
#'
#' @param sigma The standard deviation of the peak.
#'
#' @param tail A value describing any tailing desired.  If NA, no tailing is
#' applied.
#'
#' @return A vector of y values corresponding to the x values supplied.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @seealso \code{\link{lorentzCurve}}, \code{\link{makeSpec}} which uses this
#' function to make either spectra or chromatograms.
#'
#' @keywords utilities
#' @export
#' @examples
#' 
#' ### A pure Gaussian curve
#'
#' myx <- seq(0, 100, length.out = 1000) # use lots of point for resolution
#' myy <- gaussCurve(x = myx, area = 1, mu = 40, sigma = 1.5, tail = NA)
#' plot(myx, myy, type = "l", main = "Pure Gaussian Curve")
#'
#' ### Now with tailing
#' 
#' myy2 <- gaussCurve(x = myx, area = 1, mu = 40, sigma = 1.5, tail = 0.1)
#' plot(myx, myy2, type = "l", main = "Gaussian Curve with Tailing")
#' 
gaussCurve <-
function(x, area, mu, sigma, tail) {
	
	# Function to generate Gaussian curves
	# Bryan Hanson, DePauw Univ, July 2010
	# This version handles tailing

	m <- mu
	if (is.na(tail)) s <- sigma
	if (!is.na(tail)) s <- sigma*tail*x
	numerator <- exp(-1.0 * ((x - m)^2)/(2*s^2))
	denominator <- s*sqrt(2*pi)
	y <- area*numerator/denominator
	}

