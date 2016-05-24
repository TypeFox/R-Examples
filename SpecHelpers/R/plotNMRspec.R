#' Create and Plot an NMR Spectrum
#' 
#' This function simulates simple NMR spectra.  Only 1st order coupling can be
#' handled -- there is currently no capacity for doublet of doublets and
#' other such peaks.  The field strength of the "instrument" is taken into
#' account.
#' 
#' 
#' @param peaks A data frame with the following columns: delta, mult
#' (multiplicity), J, area, pw.  Multiplicity should be given by a number, so
#' use 2 for a doublet.  J is in Hz (use 0 for singlets).  pw is the peak width
#' at half-height in Hz.
#'
#' @param x.range A numeric vector of length 2 giving the ppm range desired.
#'
#' @param MHz Integer.  The operating frequency of the instrument, in MHz.
#'
#' @param ppHz Integer, but numeric works too!
#' Points per Hz: The number of data points per Hz to use in
#' calculating the spectrum (passed as argument \code{dd} to \code{makeSpec}).
#' The default (1) works well for 1H NMR spectra.  For 13C NMR spectra, where
#' the peaks are very narrow, one may need to increase the data density so that
#' enough points define the peaks (a value of 4 is a good starting point).
#' See Details.
#'
#' @param nuclei Character.  One of \code{c("1H", "13C")}. Controls the spacing
#' of the tick marks and labeling of the peaks.
#'
#' @param pkLabs Logical.  If \code{TRUE}, and \code{nuclei = 1H}, the integral
#' is drawn next to the peak.  If \code{FALSE}, no labels are drawn.
#'
#' @param lab.pos A vector of label positions as along as the number of rows in
#' \code{peaks} (the number of peaks in the spectrum).  A numeric vector
#' where 2 = left and 4 = right.  This adjusts the positions of the labels
#' to be either left or right of the peak as a way to avoid overlaps.  The
#' order must correspond to the order in \code{peaks}.
#'
#' @param plot Logical: Shall a plot be made?
#'
#' @param \dots Other parameters to be passed downstream.  These may affect
#' the plot.  You can also include \code{noise = some number} to add noise
#' (passed through to \code{makeSpec}).  In this case, warnings are raised
#' from the plotting routines, but they can be ignored.
#'
#' @return Returns a data frame of the type produced by \code{\link{makeSpec}}.
#' See there for details.  x values are in Hz.
#'
#'
#' @section Details:
#' Note that this function uses Hz internally so that the \code{x.range}, which
#' is in ppm, is multiplied by \code{Mhz} before being sent to
#' \code{\link{makeSpec}}, and once there, \code{makeSpec} will multiply it by
#' \code{ppHz}.  Thus the total data points used is \code{floor(ppHz * Mhz *
#' abs(diff(x.range)))}.  This approach ensures that peaks are not distorted
#' when changing \code{x.range} for the same \code{peak.list}.
#'
#' Note that \code{ppHz} can be numeric as well, due to the use of \code{floor}.
#' This can be useful: if you wanted your simulated NMR spectrum to be composed
#' of exactly 16384 data points as real data might be, you can call the function
#' with \code{ppHz} specified like \code{ppHz = 2^14/(12*500)} and it works! 
#'
#' @author Bryan A. Hanson, DePauw University. \email{hanson@@depauw.edu}
#'
#' @seealso \code{\link{lorentzCurve}}, \code{\link{makeSpec}}
#'
#' @keywords utilities
#' @export
#' @importFrom graphics plot axis text
#' @examples
#' 
#' ### A simulated 1H NMR spectrum
#'
#' peaks1 <- data.frame(
#' 	delta = c(1.3, 3.75, 3.9, 10.2),
#' 	mult = c(3, 4, 2, 1),
#' 	J = c(14, 14, 14, 0),
#' 	area = c(3, 2, 1, 1),
#' 	pw = c(2, 2, 2, 10))
#' 
#' res <- plotNMRspec(peaks1, x.range = c(12, 0), MHz = 500,
#' 	main = "500 MHz Simulated 1H NMR Spectrum")
#'
#' ### Compare to the same data at 200 MHz and plot together
#'
#' par(mfrow = c(2,1))
#' res <- plotNMRspec(peaks1, x.range = c(12, 0), MHz = 500,
#' 	main = "500 MHz Simulated 1H NMR Spectrum")
#' res <- plotNMRspec(peaks1, x.range = c(12, 0), MHz = 200,
#' 	main = "200 MHz Simulated 1H NMR Spectrum")
#' par(mfrow = c(1,1))
#'
#' ### Zoom in to show off
#'
#' par(mfrow = c(2,1))
#' res <- plotNMRspec(peaks1, x.range = c(4.5, 1), MHz = 500,
#' 	main = "500 MHz Simulated 1H NMR Spectrum")
#' res <- plotNMRspec(peaks1, x.range = c(4.5, 1), MHz = 200,
#' 	main = "200 MHz Simulated 1H NMR Spectrum")
#' par(mfrow = c(1,1))
#' 
#' ### A simulated 13C NMR spectrum
#'
#' # This is substantially slower due to the large
#' # chemical shift range
#' 
#' peaks2 <- data.frame(
#' 	delta = c(160, 155, 145, 143, 135, 60, 32),
#' 	mult = rep(1, 7),
#' 	J = rep(1, 7),
#' 	area = c(0.1, 0.3, 0.3, 1, 1, 0.5, 0.5),
#' 	pw = rep(1, 7))
#' 
#' res <- plotNMRspec(peaks2, x.range = c(180, 0), MHz = 200,
#' 	main = "200 MHz Simulated 13C NMR Spectrum", ppHz = 4,
#' 	pkLabs = FALSE, nuclei = "13C")
#' 
#' # Try repeating the above with ppHz = 1; note the peaks heights are not quite right
#' # as there are not enough data points to define the peak properly.
#' 
plotNMRspec <-
function(peaks, x.range = c(12, 0), MHz = 300, ppHz = 1,
	nuclei = "1H", pkLabs = TRUE, lab.pos = NULL,
	plot = TRUE, ...) {
	
	# Function to draw a theoretical 1H NMR spectrum
	# Simple, 1st order coupling only
	# Bryan Hanson, DePauw Univ, Feb 2011

	# peaks must have columns named delta, mult, J, area, pw
	# J and pw should be expressed in Hz
	# pw is peak width at half height
	# Uses makeSpec with Lorentz curves to draw the spectrum
	
	p <- peaks
	area <- ppm <- pw <- gr <- c()
	
	# Set up data to pass to makeSpec
	
	for (n in 1:nrow(p)) {
		# Get binomial coefficients for each part of multiplet
		z <- p$mult[n] - 1
		coef <- choose(z, 0:z) # a strange name for binomial!
		
		# Calc areas
		a <- coef * p$area[n]/sum(coef)
		
		# Calc chemical shifts
		no.pks <- length(coef)
		m <- jSeq(no.pks) 
		m <- p$J[n]*m
		m <- m + p$delta[n]*MHz # positions in Hz
		gr <- c(gr, rep(n, no.pks))
		area <- c(area, a)
		ppm <- c(ppm, m)
		pw <- c(pw, rep(p$pw[n], length(m)))
		}

	ans <- data.frame(area = area, x0 = ppm, gamma = pw/2)

	# p.ppm = points per ppm This way the data density
	# is independent of the requested x.range
	# This is esp. important for 13C since the peaks are narrow
	
	spec <- makeSpec(peak.list = ans, plot = FALSE,
		x.range = x.range*MHz, type = "lorentz", dd = ppHz, ...)


	# Labels are actually the hardest part!
	
	if (pkLabs) { # Figure out where the label will go...
		
		spec2 <- spec[-c(1,2),]
		grps <- as.factor(gr)

		y.pos <- NA_real_
		for (i in 1:length(levels(grps))) {
			w <- which(grps == levels(grps)[i])
			y.pos <- c(y.pos, max(spec2[w,]))
			}
		y.pos <- y.pos[-1]
#		cat("y.pos original =", y.pos, "\n")
		
		# Next part adds an extra offset if labels are close together
		# A better approach than this would be to calc exactly where
		# the label is and move high enough
		extra <- 0.5*min(y.pos)
		y.off <- rep(0, (length(y.pos)-1))
#		cat("y.off initially =", y.off, "\n")
		gap <- abs(diff(peaks$delta))
#		cat("gap =", gap, "\n")
		if (any(gap <= 0.3)) {
			for (n in 1:length(gap)) {
				if (gap[n] <= 0.3) y.off[n] <- 1
				}	
			}
#		cat("y.off final =", y.off, "\n")
		y.off <- c(0, extra*y.off)
		y.pos <- y.pos + y.off
#		cat("y.pos final =", y.pos, "\n")
		# End of extra offset calculation
		
		if (nuclei == "1H") labs <- paste(p$area, "H", sep = " ")
		if (nuclei == "13C") labs <- as.character(p$delta)
		}

	# Now do the plotting
	
	if (plot) {
		if (!pkLabs) y.lim <- range(spec[-c(1,2),])
		if (pkLabs) y.lim <- c(0, max(y.pos)*1.1) # Increase ylim a bit so integral labels are not cut off
			
		plot(x = spec[1,]/MHz, y = spec[2,],
			type = "l", xlim = x.range, xlab = "chemical shift, ppm",
			ylab = "", axes = FALSE, ylim = y.lim, ...)
			
		if (nuclei == "1H") {
			# lbl <- rep(NA_real_, 76)
			# lbl[which(1:76 %% 5 == 1L)] <- round(0:15, 1)
			axis(side = 1, at = seq(0, 15, by = 0.2), tcl = -0.25, labels = FALSE)
			axis(side = 1, at = 0:15)
			}
			
		if (nuclei == "13C") {
			axis(side = 1, at = seq(0, 200, by = 10), tcl = -0.25, labels = FALSE)
			axis(side = 1, at = seq(0, 200, by = 20))
			}
		if (pkLabs) { # Now that the plotting window is open add the labels
			lp = 2
			if (!is.null(lab.pos)) lp <- lab.pos # need to add translation of L, R to 2, 4
			text(x = p$delta, y = y.pos, labels = labs, col = "red", pos = lp, offset = 0.5)
			}
	}
		
	return(spec)
	}
