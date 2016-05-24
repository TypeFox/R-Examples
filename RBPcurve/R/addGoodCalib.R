#' @title Visualizes a measure for good calibration on the RBP curve.
#'
#' @description The integral of the RBP curve is a measure for good calibration.
#' If the sum of the two integrals (below and above the RBP curve) is close to 0,
#' good calibration is satisfied and the prevalence is close to the average predicted probabilities.
#'
#' @template arg_obj
#' @template arg_plotvalues
#' @template arg_showinfo
#' @param col [\code{vector(1)}]\cr
#'   Color for filling the polygon, as in \code{\link{polygon}}.
#'   Default is \dQuote{grey}.
#' @param border [\code{vector(1)}]\cr
#'   Color to draw the borders, as in \code{\link{polygon}}.
#'   Default is \code{NA} to omit borders.
#' @param ... [any]\cr
#'   Passed to \code{\link{polygon}}.
#' @template ret_invnull
#' @export
addGoodCalib = function(obj, plot.values = TRUE, show.info = TRUE,
  col = grDevices::rgb(0, 0, 0, 0.25), border = NA, ...) {

  # Check arguments
  assertClass(obj, "RBPObj")
  assertFlag(plot.values)
  assertVector(col, len = 1L)
  assertVector(border, len = 1L)

  # Store values of obj
  x1 = obj$axis.x
  y1 = obj$axis.y
  eps = obj$y - obj$pred

  # Highlights the integral below the RBP curve
  polygon(c(min(x1), x1, max(x1)), c(0, y1, 0), col = col, border = border, ...)

  # Integral below and above the RBP curve
  below = round(sum(eps[obj$y == 0]) / obj$n, 4L)
  above = round(sum(eps[obj$y == 1]) / obj$n, 4L)
  
  # Add the values of the integral below and above the RBP curve into the current plot?
  if (plot.values) {
    text(x1[sum(obj$y == 0)], 0, adj = 1:0, labels = below)
    text(x1[sum(obj$y == 0)+1], 0, adj = 0:1, labels = above)
  }

  # Show message
  if (show.info) {
    messagef("Integral below the RBP curve: %s", below)
    messagef("Integral above the RBP curve: %s", above)
  }

  return(invisible(NULL))
}

