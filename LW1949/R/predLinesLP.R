#' Add Litchfield and Wilcoxon Predictions to a Plot
#'
#' Add predictions from a Litchfield and Wilcoxon model fit to
#' a plot of the results of a dose-effect experiment on the log10-probit scale.
#' @param fit
#'   A list of length three containing the result of a Litchfield and
#'     Wilcoxon model fit, typically the output from \code{\link{LWestimate}}.
#' @param ...
#'   Additional arguments to \code{\link{abline}} and \code{\link{lines}}.
#' @return
#'   A solid fitted line is added to the plot.  Dashed lines are added to the
#'     plot representing the \strong{horizontal} 95% confidence intervals
#'     for the predicted dose to elicit a given percent affected.
#' @export
#' @import
#'   graphics
#' @seealso  
#'   \code{\link{plotDELP}}, \code{\link{plotDE}}, \code{\link{predLines}}
#' @examples
#' dose <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' ntested <- rep(8, 5)
#' nalive <- c(1, 4, 4, 7, 8)
#' mydat <- dataprep(dose=dose, ntot=ntested, nfx=nalive)
#' plotDELP(mydat)
#' myfit <- LWestimate(fitLWauto(mydat), mydat)
#' predLinesLP(myfit)

predLinesLP <- function(fit, ...) {
  ys <- c(seq(0.1, 0.9, 0.1), 1:99, seq(99.1, 99.9, 0.1))
  lc <- predlinear(ys, fit)
  abline(fit$params, lwd=2, ...)
  lines(log10(lc[, "upper"]), probit(lc[, "pct"]/100), lty=2, ...)
  lines(log10(lc[, "lower"]), probit(lc[, "pct"]/100), lty=2, ...)
}
