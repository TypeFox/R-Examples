#' Plot Dose-Effect Experiments
#'
#' Plot dose-effect experiments on the log10-probit scale.
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
#' @param ...
#'   Additional arguments to \code{\link{plot}}.
#' @export
#' @import
#'   graphics
#' @seealso  
#'   \code{\link{predLinesLP}}, \code{\link{plotDE}}, \code{\link{predLines}}
#' @examples
#' dose <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' ntested <- rep(8, 5)
#' nalive <- c(1, 4, 4, 7, 8)
#' mydat <- dataprep(dose=dose, ntot=ntested, nfx=nalive)
#' plotDELP(mydat)

plotDELP <- function(DEdata, xlab="Dose", ylab="Affected  (%)",
  xlim=range(DEdata$dose[DEdata$dose>0], na.rm=TRUE), ylim=c(0.1, 99.9), ...) {
  if (!is.data.frame(DEdata)) stop("DEdata must be a data frame.")
  if (any(is.na(match(c("dose", "pfx", "log10dose", "bitpfx", "fxcateg"),
    names(DEdata))))) {
    stop("DEdata must include at least five variables:",
      "dose, pfx, log10dose, bitpfx, fxcateg.")
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
  if (any(ylim <= 0) | any(ylim >= 100)) stop("ylim must be between 0 and 100.")

  xtix <- prettylog(c(DEdata$dose, xlim), 1:9, 5)
  xlim <- range(log10(c(DEdata$dose[DEdata$dose>0], xlim)), na.rm=TRUE)

  plot(DEdata$log10dose, DEdata$bitpfx, type="n", xlim=xlim,
    ylim=probit(ylim/100), axes=F, xlab=xlab, ylab=ylab, ...)

  # background grid and axes
  abline(v=log10(xtix), lwd=2, col="lightgray")
  axis(1, at=log10(xtix), labels=xtix)

  ytix1 <- c(seq(0.1, 0.9, 0.1), seq(1, 9, 1), seq(10, 90, 10), seq(91, 99, 1),
    seq(99.1, 99.9, 0.1))
  ytix2 <- c(0.1, 1, 10, 50, 90, 99, 99.9)
  abline(h=probit(ytix1/100), col="lightgray")
  abline(h=probit(ytix2/100), lwd=2, col="lightgray")
  axis(2, at=probit(ytix2/100), labels=ytix2, las=1)
  box()

  # observed points
  points(DEdata$log10dose, probit(constrain(DEdata$pfx, ylim/100)), pch=16,
    cex=1.5)
  points(DEdata$log10dose[DEdata$fxcateg!=50],
    probit(constrain(DEdata$pfx[DEdata$fxcateg!=50], ylim/100)),
    pch=16, col="white")
  }
