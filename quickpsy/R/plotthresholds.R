#' Plot the thresholds
#'
#' \code{plotthresholds} plot the thresholds.
#' @param qp output from quickpsy.
#' @param x Name of the variable to displayed in the x-axis.
#' @param panel Name of the variable to be split in panels.
#' @param xpanel Name of the variable to be split in horizontal panels.
#' @param ypanel Name of the variable to be split in vertical panels.
#' @param color Name of the variable codded by color.
#' @param geom If \code{'bar'} displays bars.
#' @param sizeerrorbar Line width of the error bars.
#' If \code{'point'} displays points (default is 'bar').
#' @param ci If \code{FALSE} confidence intervals are not plotted
#' (default is \code{TRUE}).
#' @seealso  \code{\link{plotthresholds_}}
#' @examples
#' library(MPDiR) # contains the Vernier data
#' fit <- quickpsy(Vernier, Phaseshift, NumUpward, N,
#'                 grouping = .(Direction, WaveForm, TempFreq), B = 10)
#' plotthresholds(fit)
#' plotthresholds(fit, x = WaveForm)
#' plotthresholds(fit, xpanel = Direction)
#' plotthresholds(fit, color = Direction, ypanel = WaveForm, geom = 'point')
#' @export
plotthresholds <- function(qp, x = NULL, panel = NULL, xpanel = NULL,
                           ypanel = NULL, color = NULL, geom = 'bar', ci = T,
                           sizeerrorbar = 1) {

  if (!missing(x)) x <- deparse(substitute(x))
  if (!missing(panel)) panel <- deparse(substitute(panel))
  if (!missing(xpanel)) xpanel <- deparse(substitute(xpanel))
  if (!missing(ypanel)) ypanel <- deparse(substitute(ypanel))
  if (!missing(color)) color <- deparse(substitute(color))

 plotthresholds_(qp, x, panel, xpanel, ypanel, color, geom, ci,sizeerrorbar)
}
