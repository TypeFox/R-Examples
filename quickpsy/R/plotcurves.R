#' Plot the curves
#'
#' \code{plotcurves} plot the curves.
#' @param qp output from quickpsy
#' @param panel Name of the variable to be split in panels.
#' @param xpanel Name of the variable to be split in horizontal panels.
#' @param ypanel Name of the variable to be split in vertical panels.
#' @param color Name of the variable codded by color.
#' @param averages If \code{FALSE} averaged probabilities are not plotted
#' (default is \code{TRUE}).
#' @param curves If \code{FALSE} curves are not plotted
#' (default is \code{TRUE})
#' @param thresholds If \code{FALSE} thresholds  are not plotted
#' (default is \code{TRUE})
#' @param ci If \code{FALSE} confidence intervals are not plotted
#' (default is \code{TRUE})
#' @seealso \code{\link{plotcurves_}}
#' @examples
#' library(MPDiR) # contains the Vernier data
#' fit <- quickpsy(Vernier, Phaseshift, NumUpward, N,
#'                 grouping = .(Direction, WaveForm, TempFreq), B = 10)
#' plotcurves(fit)
#' plotcurves(fit, xpanel = Direction)
#' plotcurves(fit, xpanel = Direction, color = WaveForm, ci = FALSE)
#' @export
plotcurves <- function(qp, panel = NULL, xpanel = NULL, ypanel = NULL,
                       color = NULL, averages = T, curves = T, thresholds = T,
                       ci = T) {

  if (!missing(panel)) panel <- deparse(substitute(panel))
  if (!missing(xpanel)) xpanel <- deparse(substitute(xpanel))
  if (!missing(ypanel)) ypanel <- deparse(substitute(ypanel))
  if (!missing(color)) color <- deparse(substitute(color))

  plotcurves_(qp, panel, xpanel, ypanel, color, averages, curves,
              thresholds, ci)
}
