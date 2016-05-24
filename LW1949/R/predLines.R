#' Add Litchfield and Wilcoxon Predictions to a Plot
#'
#' Add predictions from a Litchfield and Wilcoxon model fit to
#' a plot of the results of a dose-effect experiment on the arithmetic scale.
#' @param fit
#'   A list of length three containing the result of a Litchfield and
#'     Wilcoxon model fit, typically the output from \code{\link{LWestimate}}.
#' @param ...
#'   Additional arguments to \code{\link{lines}}.
#' @return
#'   A solid fitted line is added to the plot.  Dashed lines are added to the
#'     plot representing the \strong{horizontal} 95% confidence intervals
#'     for the predicted dose to elicit a given percent affected.
#' @export
#' @import
#'   graphics
#' @seealso  
#'   \code{\link{plotDE}}, \code{\link{plotDELP}}, \code{\link{predLinesLP}}
#' @examples
#' dose <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' ntested <- rep(8, 5)
#' nalive <- c(1, 4, 4, 7, 8)
#' mydat <- dataprep(dose=dose, ntot=ntested, nfx=nalive)
#' plotDE(mydat)
#' myfit <- LWestimate(fitLWauto(mydat), mydat)
#' predLines(myfit)

predLines <- function(fit, ...) {
  ys <- c(seq(0.1, 0.9, 0.1), 1:99, seq(99.1, 99.9, 0.1))
  lc <- predlinear(ys, fit)
  lines(lc[, "ED"], lc[, "pct"], lwd=2, ...)
  lines(lc[, "upper"], lc[, "pct"], lty=2, ...)
  lines(lc[, "lower"], lc[, "pct"], lty=2, ...)
}
