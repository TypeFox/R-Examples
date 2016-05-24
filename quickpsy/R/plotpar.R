#' Plot the values of the parameters
#'
#' \code{plotpar} plot the values of the parameters.
#' @param qp output from quickpsy.
#' @param x Name of the variable to displayed in the x-axis.
#' @param panel Name of the variable to be split in panels.
#' @param xpanel Name of the variable to be split in horizontal panels.
#' @param ypanel Name of the variable to be split in vertical panels.
#' @param color Name of the variable codded by color.
#' @param geom If \code{'bar'} displays bars.
#' If \code{'point'} displays points (default is \code{'bar'}).
#' @param ci If \code{FALSE} confidence intervals are not plotted
#' (default is \code{TRUE}).
#' @seealso  \code{\link{plotpar_}}
#' @examples
#' library(MPDiR) # contains the Vernier data
#' fit <- quickpsy(Vernier, Phaseshift, NumUpward, N,
#'                 grouping = .(Direction, WaveForm, TempFreq), B = 10)
#' plotpar(fit)
#' plotpar(fit, x = WaveForm)
#' plotpar(fit, xpanel = Direction)
#' plotpar(fit, color = Direction)
#' plotpar(fit, color = Direction, ypanel = WaveForm, geom = 'point')
#' @export
plotpar <- function(qp, x = NULL, panel = NULL, xpanel = NULL,
                           ypanel = NULL, color = NULL, geom = 'bar', ci = T) {

  if (!missing(x)) x <- deparse(substitute(x))
  if (!missing(panel)) panel <- deparse(substitute(panel))
  if (!missing(xpanel)) xpanel <- deparse(substitute(xpanel))
  if (!missing(ypanel)) ypanel <- deparse(substitute(ypanel))
  if (!missing(color)) color <- deparse(substitute(color))

 plotpar_(qp, x, panel, xpanel, ypanel, color, geom)
}
