#' Draw a Chromatogram or Spectrum
#' 
#' This function creates a chromatogram or spectrum from a list of appropriate
#' parameters describing the peaks.  The individual curves are computed using
#' the mathematical definition of either a Gaussian curve, possibly with
#' tailing, or a Lorentzian curve.  Gaussian curves are appropriate for
#' simulating chromatograms or UV-Vis spectra, while Lorentzians are
#' used for simulating NMR peaks.  The function computes the individual curves
#' as well as their sum (which is the whole chromatogram or spectrum).  A plot
#' can be made, which may display the separate underlying curves.  If you want
#' to draw NMR spectra, use \code{\link{plotNMRspec}} which is a much more
#' natural interface to this function.
#' 
#' 
#' @param peak.list For a Gaussian curve, a data frame with the following
#' columns: mu, sd, area, tail.  mu is the retention time (or center
#' frequency).  sd is the standard deviation (or peak width).  area is the area
#' under the peak.  tail is the tailing parameter - use NA when a pure Gaussian
#' with no tailing is desired.  One row of the data frame contains data related
#' to one peak.
#' 
#' For a Lorentzian curve, a data frame with the following columns: x0, area,
#' gamma.  x0 is the center frequency or chemical shift.  gamma is the half the
#' peak width at half-height.  area is the area under the peak.
#'
#' @param x.range A numeric vector of length 2 giving the retention time range
#' (or frequency range) desired.  Must make sense in light of the peak list
#' given (i.e. a wider range, possibly much wider depending up the values of
#' \code{sd} and \code{tail}), as these broaden the peaks.
#'
#' @param plot Logical; if TRUE, a plot is produced.
#'
#' @param curves Logical; if TRUE, the individual curves are plotted (provided
#' \code{plot = TRUE}.  Not very useful for NMR spectra, but great for showing,
#' for instance, how shoulders arise on peaks in a chromatogram.
#'
#' @param type A character string.  Use "gauss" to generate Gaussian curves
#' (for chromatograms, or UV-Vis spectra).  Use "lorentz" to generate
#' Lorentzian curves as found in NMR spectra.
#'
#' @param noise A number giving the amount of noise to be added to the
#' individual curves (the net spectrum has the noise from the individual
#' spectra, it has no additional noise added to it).  Value corresponds to the
#' argument \code{factor} in function \code{jitter}.
#'
#' @param dd The density of data points per unit of \code{x.range}.  The total
#' number of data points used to create the spectrum or chromatogram is
#' \code{dd*abs(diff(x.range))} and thus it also depends on the units of
#' \code{x.range}.  This approach ensures that peaks are not distorted when
#' changing \code{x.range} for the same \code{peak.list}.
#'
#' @param \dots Additional arguments to be passed downstream.
#'
#' @return A matrix containing the x values (retention times or
#' frequencies) in the first row, and the complete chromatogram (spectrum) in
#' the second row.  Additional rows contain chromatograms (spectra) of the
#' individual components.  The row names of the data frame are character
#' strings describing the chromatogram (spectrum) in that row.  The matrix
#' contains \code{dd*abs(diff(x.range))} columns.
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @seealso \code{\link{gaussCurve}}, \code{\link{lorentzCurve}},
#' \code{\link{plotNMRspec}} and \code{\link{plot2DNMRspec}}, the preferred
#' interfaces for drawing NMR spectra.
#' @keywords utilities
#' @export
#' @importFrom graphics plot lines
#' @examples
#' 
#' ### A simple chromatogram
#' 
#' chrom <- data.frame(mu = c(2, 5, 11), sd = c(0.5, 1, 2),
#' area = c(1, 0.5, 1), tail =  c(NA, NA, 0.1))
#' ex1 <- makeSpec(chrom, x.range = c(0, 20), plot = TRUE, curves = TRUE,
#' dd = 5, main = "Chromatogram with Underlying Pure Curves")
#'
#' ### Faux ethyl group NMR with J = 0.1 ppm.
#' # Note that a much better
#' # NMR spectrum can be generated using plotNMRspec which also uses
#' # a more natural input format
#' #
#' spec <- data.frame(mu = c(3.5, 3.4, 3.3, 3.2, 1.4, 1.3, 1.2),
#' sd = rep(0.01, 7), tail =  rep(NA, 7),
#' area = c(1, 3, 3, 1, 1, 2, 1) * c(0.5, 0.5, 0.5, 0.5, 0.66, 0.66, 0.66))
#' ex2 <- makeSpec(spec, x.range = c(5, 0), plot = TRUE, curves = FALSE,
#' dd = 100, main = "Simulated 1H NMR of an Ethyl Group")
#' 
makeSpec <-
function(peak.list, x.range, plot = TRUE, curves = FALSE,
	type = "gauss", noise = 0, dd = 1, ...) {

# Function to generate sample spectra or chromatograms
# Bryan Hanson, DePauw Univ, July 2010

	# For Gaussian curves:
	# peak.list must contain area, mu, sd, tail
	# tailing is handled by making sd a function of x (time)
	# tailing can be disabled entirely with tail = NA

	# For Lorentzian curves:
	# peak.list must contain area, x0, gamma
	# tailing is not relevant
	# Conceptually, x0 ~ mu, gamma ~ sd

	# dd = data density per x.range unit
	# ndp = no. data points (total)
	# Note that plotNMRspec actually passes x.range in Hz not ppm

	ndp <- floor(dd*abs(diff(x.range)))
#	cat("makeSpec says: ndp = ", ndp, "\n")

	if (type == "gauss") {		
		pl <- peak.list
		ns <- length(pl$mu) # ns = no. of spec
		if (is.null(pl$tail)) pl$tail <- rep(NA, ns)
	
		# create x-data, initialize y-data
		# y.mat will hold each spectrum separately

		x <- seq(from = x.range[1], to = x.range[2], length.out = ndp)
		y.mat <- matrix(data = NA_real_, nrow = ns, ncol = ndp)
	
		for (n in 1:ns) {
			y.mat[n,] <- gaussCurve(x = x, area = pl$area[n], mu = pl$mu[n],
				sigma = pl$sd[n], tail = pl$tail[n])
			}

		rn <- list()
		for (n in 1:ns) {
			rn[n] <- paste("area", pl$area[n], "mu", pl$mu[n], "sigma", pl$sd[n], "tail", pl$tail[n], sep = " ")
			}
		}

	if (type == "lorentz") {		
		pl <- peak.list
		ns <- length(pl$x0) # ns = no. of spec
		x <- seq(from = x.range[1], to = x.range[2], length.out = ndp)
		y.mat <- matrix(data = NA_real_, nrow = ns, ncol = ndp)
				
		for (n in 1:ns) {
			y.mat[n,] <- lorentzCurve(x = x, area = pl$area[n],
				x0 = pl$x0[n], gamma = pl$gamma[n])
			}

		rn <- list()
		for (n in 1:ns) {
			rn[n] <- paste("area", pl$area[n], "x0", pl$x0[n], "gamma", pl$gamma[n], sep = " ")
			}
		}

	dimnames(y.mat)[[1]] <- rn
	
	if (!noise == 0) { y.mat <- jitter(y.mat, factor = noise)}

	y.sum <- colSums(y.mat)
	all <- rbind(x, y.sum, y.mat)

	if (plot) {
		plot(x, y.sum, type = "l", lwd = 2, col = "black", xlim = x.range, ...)
		if (curves) for (n in 1:ns) lines(x, y.mat[n,], lwd = 1.0, col = "blue")
		}
	
	return(all)
	}

