#' Draws a caterpillar plot (in MLwiN style).
#'
#' A convenient wrapper for the \code{\link[graphics]{plot}} function with the addition
#' of error bars, e.g. to create caterpillar plots.
#'
#' @param y A numerical vector specifying the \code{y} coordinates
#' (e.g. the means or medians); see \code{\link[graphics]{plot.default}}.
#' @param x A numerical vector specifying the \code{x} coordinates; see \code{\link[graphics]{plot.default}}.
#' @param qtlow A numerical vector (e.g. of lower-quantiles) to be used to plot lower error
#' bars.
#' @param qtup A numerical vector (e.g. of upper-quantiles) to be used to upper plot error
#' bars.
#' @param xlab A label for the \code{x} axis. This is empty by default; see \code{\link[graphics]{plot.default}}.
#' @param ylab A label for the \code{y} axis. This is empty by default; see \code{\link[graphics]{plot.default}}.
#' @param xlim The \code{x} limits \code{(x1, x2)} of the plot. Note that
#' \code{x1 > x2} is allowed and leads to a `reversed axis'. The default value,
#' \code{NULL}, indicates that the range of the finite values to be plotted
#' should be used; see \code{\link[graphics]{plot.default}}.
#' @param ylim The y limits of the plot; see \code{\link[graphics]{plot.default}}.
#' @param main A main title for the plot; see \code{\link[graphics]{plot.default}}.
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol, U.K.
#'
#' @seealso \code{\link{caterpillarR}}
#'
#' @examples
#'
#' \dontrun{
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved where R2MLwiN defaults to:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/')
#'
#' # Example using tutorial dataset
#' data(tutorial, package = 'R2MLwiN')
#' (mymodel <- runMLwiN(normexam ~ 1 + (1 | school) + (1 | student),
#'                      estoptions = list(resi.store = TRUE),
#'                      data = tutorial))
#'
#' # For each school, calculate the CIs...
#' residuals <- mymodel@@residual$lev_2_resi_est_Intercept
#' residualsCI <- 1.96 * sqrt(mymodel@@residual$lev_2_resi_var_Intercept)
#' residualsRank <- rank(residuals)
#' rankno <- order(residualsRank)
#'
#' caterpillar(y = residuals[rankno], x = 1:65, qtlow = (residuals - residualsCI)[rankno],
#'            qtup = (residuals + residualsCI)[rankno], xlab = 'Rank', ylab = 'Intercept')
#' }
#'
#' @export
caterpillar <- function(y, x, qtlow, qtup, xlab = "", ylab = "", xlim = NULL, ylim = NULL, main = "") {
  # This function draws a caterpillar plot
  if (is.null(ylim)){
    ylim = c(min(qtlow), max(qtup))
  }
  plot(x, y, xlim = xlim, ylim = ylim, pch = 15, xlab = xlab, ylab = ylab, main = main)
  points(x, qtlow, pch = 24, bg = "grey")
  points(x, qtup, pch = 25, bg = "grey")
  for (i in 1:length(x)) {
    lines(rep(x[i], 2), c(qtlow[i], qtup[i]))
  }
}
