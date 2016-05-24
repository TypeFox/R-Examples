#' Function to estimate the approximate local minima and maxima of melting
#' curve data.
#' 
#' The \code{mcaPeaks()} is used to estimate the approximate local minima and
#' maxima of melting curve data. This can be useful to define a temperature
#' range for melting curve analysis, quality control of the melting curve or to
#' define a threshold of peak heights. Melting curves may consist of multiple
#' significant and insignificant melting peaks. \code{mcaPeaks()} uses
#' estimated the temperatures and fluorescence values of the local minima and
#' maxima. The original data remain unchanged and only the approximate local
#' minima and maxima are returned. \code{mcaPeaks()} takes modified code
#' proposed earlier by Brian Ripley
#' (https://stat.ethz.ch/pipermail/r-help/2002-May/021934.html).
#' 
#' 
#' @param x \code{x} is the column of a data frame for the temperature.
#' @param y \code{y} is the column of a data frame for the fluorescence values.
#' @param span \code{span} is used to adjust the window span.
#' @return \item{p.min }{returns a \code{data.frame} with the temperatures ("T
#' (minima)") and the fluorescence ("F (minima)") for the local minima of the
#' processed data. }
#' 
#' \item{p.max }{returns a \code{data.frame} with the temperatures ("T
#' (minima)") and the fluorescence ("F (minima)") for the local maxima of the
#' processed data. }
#' @author Stefan Roediger
#' @seealso \code{\link{mcaSmoother}}, \code{\link{smooth.spline}}
#' @keywords smooth peaks
#' @examples
#' 
#' # First Example
#' data(DMP)
#' # Smooth and Min-Max normalize melting curve data with mcaSmoother()
#' tmp <- mcaSmoother(DMP[, 1], DMP[,6], minmax = TRUE, n = 2)
#' 
#' # Extract the first derivative melting curve data
#' tmp.d <- diffQ(tmp, verbose = TRUE)$xy
#' 
#' # Determine the approximate local minima and maxima of a curve
#' peak.val <-mcaPeaks(tmp.d[, 1], tmp.d[, 2])
#' peak.val
#' 
#' # Plot the first derivative melting curve
#' # Add a horizontal line and points for the approximate local minima 
#' # and maxima to the plot
#' plot(tmp.d, type = "l", 
#'      main = "Show the approximate local minima and maxima of a curve")
#'   abline(h = 0)
#'   points(peak.val$p.min, col = 1, pch = 19)
#'   points(peak.val$p.max, col = 2, pch = 19)
#'   legend(25, 0.08, c("Minima", "Maxima"), col = c(1,2), pch = c(19,19))
#' 
#' # Second example
#' # Signifcant peaks can be distinguished by peak hight
#' plot(tmp.d, type = "l", 
#'       main = "Show the approximate local minima and maxima of a curve")
#'   abline(h = 0)
#'   points(peak.val$p.min, col = 1, pch = 19)
#'   points(peak.val$p.max, col = 2, pch = 19)
#'   legend(25, 0.08, c("Minima", "Maxima"), col = c(1,2), pch = c(19,19))
#' 
#' # Test which local maxima peak is above the median + 3 * Median Absolute 
#' # Add a threshold (th) line to the plot
#' th.max <- median(peak.val$p.max[, 2]) + 3 * mad(peak.val$p.max[, 2])
#' abline(h = th.max, col = 3)
#' 
#' # add the approximate temperatures of the local minima and 
#' # maxima to the plot
#' T.val <- c(peak.val$p.min[, 1], peak.val$p.max[, 1])
#' text(T.val, rep(-0.05, length(T.val)), round(T.val,0))
#' 
#' # Use a temperature range from the plot to calculate the Tm of 
#' # the maximum Trange is used between 37 and 74 degree Celsius
#' 
#' tmp <- mcaSmoother(DMP[, 1], DMP[, 6], minmax = TRUE, Trange = c(37,74), 
#' 		    n = 2)
#' # Tm 48.23, fluoTm 0.137
#' diffQ(tmp, fct = max, plot = TRUE)
#' 
#' 
#' @export mcaPeaks
mcaPeaks <- function(x, y, span = 3) {
  # Take the x and y values from the object of type data.frame.
  
  old.warn <- options("warn")[["warn"]]
  options(warn = -1)
  # Test if x and y exist.
  if (is.null(x)) 
      stop("Enter temperature")
  if (is.null(y)) 
      stop("Enter fluorescence data")
  # Test if x and y have the same length.
  if (length(x) != length(y)) 
      stop("Use temperature and fluorescence data with same number of elements")
  # Test if span is odd.
  if (span%/%2 != 1) {
      message("span must be odd. Automatically set to span = 3")
      span <- 3
  }

  # smooth the input data slightly with a spline
  input.y <- smooth.spline(x, y)[["y"]]

  # Estimate the local minima and local maxima
  # autor Brian Ripley, https://stat.ethz.ch/pipermail/r-help/2002-May/021934.html
  # modified for the MBmca package to find the approximate maxima and minima of
  # the melting peak data
  s <- span%/%2

  z.max <- embed(input.y, span)
  v.max <- max.col(z.max) == 1 + s
  maxima <- c(rep(FALSE, s), v.max)
  maxima <- maxima[1L:(length(maxima) - s)]

  z.min <- embed(input.y * -1, span)
  v.min <- max.col(z.min) == 1 + s
  minima <- c(rep(FALSE, s), v.min)
  minima <- minima[1L:(length(minima) - s)]

  out <- data.frame(x, input.y)

  p.min <- data.frame(out[which(minima == TRUE), 1], 
		      out[which(minima == TRUE), 2]
		      )
  colnames(p.min) <- c("T (minima)", "F (minima)")

  p.max <- data.frame(out[which(maxima == TRUE), 1], 
		      out[which(maxima == TRUE), 2]
		      )
  colnames(p.max) <- c("T (maxima)", "F (maxima)")
  
  #restore old warning value
  options(warn = old.warn)
  
  list(p.min = p.min, p.max = p.max)
} 
