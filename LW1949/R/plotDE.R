#' Plot Dose-Effect Experiments
#'
#' Plot dose-effect experiments on the arithmetic scale.
#' @param DEdata
#'   A data frame of dose-effect data (typically, the output from
#'     \code{\link{dataprep}}) containing at least five variables:
#'     dose, pfx, log10dose, bitpfx, fxcateg.
#' @param xlab
#'   A character scalar, the title for the dose (x) axis, default "Dose".
#' @param ylab
#'   A character scalar, the title for the affected (y) axis,
#'     default "Affected  (\%)".
#' @param xlim
#'   A numeric vector of length two giving the x coordinate range for
#'     dose, default range(DEdata$dose, na.rm=TRUE).
#' @param ylim
#'   A numeric vector of length two giving the y coordinate range for
#'     affected (\%), default c(0.1, 99.9).
#'   Observed effects beyond this range will be plotted at the limits of this
#'     range using an open symbol.
#' @param ref
#'   A numeric vector specifying horizontal reference lines to be added
#'     to the plot, default c(0, 50, 100).
#' @param ...
#'   Additional arguments to \code{\link{plot}}.
#' @export
#' @import
#'   graphics
#' @seealso  
#'   \code{\link{predLines}}, \code{\link{plotDELP}}, \code{\link{predLinesLP}}
#' @examples
#' dose <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' ntested <- rep(8, 5)
#' nalive <- c(1, 4, 4, 7, 8)
#' mydat <- dataprep(dose=dose, ntot=ntested, nfx=nalive)
#' plotDE(mydat)

plotDE <- function(DEdata, xlab="Dose", ylab="Affected  (%)",
  xlim=range(DEdata$dose, na.rm=TRUE), ylim=c(0, 100), ref=c(0, 50, 100), ...) {
  if (!is.data.frame(DEdata)) stop("DEdata must be a data frame.")
  if (any(is.na(match(c("dose", "pfx", "fxcateg"), names(DEdata))))) {
    stop("DEdata must include at least three variables:",
      "dose, pfx, fxcateg.")
  }
  if (!is.character(xlab) | length(xlab)!=1) {
    stop("xlab must be a character scalar")
  }
  if (!is.character(ylab) | length(ylab)!=1) {
    stop("ylab must be a character scalar")
  }
  if (!is.numeric(ylim) | length(ylim)!=2) {
    stop("ylim must be a numeric vector of length 2.")
  }
  if (any(ylim < 0) | any(ylim > 100)) {
    stop("ylim must be between 0 and 100, inclusive.")
  }

  y <- 100*DEdata$pfx
  plot(DEdata$dose, y, type="n", xlim=xlim,
    ylim=ylim, xlab=xlab, ylab=ylab, las=1, ...)

  # reference lines
  abline(h=ref, lwd=2, col="lightgray")

  # observed points
  points(DEdata$dose, y, pch=16, cex=1.5)
  sel <- DEdata$fxcateg!=50
  points(DEdata$dose[sel], y[sel], pch=16, col="white")
}
