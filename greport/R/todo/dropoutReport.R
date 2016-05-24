#' Dropout Report
#'
#' Generate a survival plot for subjects remaining in the study.
#'
#' @param d.dropout numeric vector. Dropout date.
#' @param dropout numeric vector. Indicator variable for dropout.
#' @param treat factor vector. Vector of treatment group for each record.
#' @param time.inc numeric. See \code{\link[rms]{survplot}}.
#' @param ylim numeric vector. See \code{\link[rms]{survplot}}.
#' @param panel character. Name for panel, defaults to \sQuote{dropout}.
#' @param what character. Name of study, defaults to \sQuote{study}.
#' @param \dots additional arguments, passed to \code{\link[rms]{survplot}}.
#' @export
#' @examples
#' \dontrun{
#'   d.d <- sample(1:10, 200, replace=TRUE, prob=c(rep(0.03,9), 0.73))
#'   dropout <- as.numeric(d.d < 10)
#'   dropoutReport(d.d, dropout, as.factor(sample(c('A','B'), 200, replace=TRUE)), time.inc=2)
#' }

dropoutReport <- function(d.dropout, dropout, treat,
                          time.inc=NULL, ylim=c(0,1), panel="dropout", what="study", ...) {
  ### function for capitalizing the first letter of each word
  ### borrowed from function "toupper" help
  .simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
  }
  S <- if(length(dropout)) Surv(d.dropout, dropout) else
      Surv(d.dropout)
  openPanel <- paste("O", panel, sep="")
  startPlot(openPanel, h=4)
  f <- survfit.formula(S ~ treat)
  d <- data.frame(treat)
  d$S <- S
  yl <- paste("Fraction Remaining in",.simpleCap(what))
  lwd <- c(1,2); lty=c(1,1); col=gray(c(0,.7))
  if(length(time.inc))
    survplot.survfit(f, time.inc=time.inc, n.risk=TRUE, conf='none', ylab=yl,
      lwd=lwd, lty=lty, col=col, ylim=ylim, ...) else
      survplot.survfit(f, conf='none', ylab=yl, lwd=lwd, lty=lty, col=col,
               ylim=ylim, n.risk=TRUE, ...)

  endPlot()
  startPlot(panel, h=4)
  f <- survfit.formula(S ~ treat, data=d)
  if(length(time.inc))
    survplot.survfit(f, time.inc=time.inc, n.risk=TRUE, conf='none', ylab=yl,
             lwd=lwd, lty=lty, col=col, label.curves=FALSE, ylim=ylim, ...)
      else survplot.survfit(f, conf='none', ylab=yl, n.risk=TRUE,
                    lwd=lwd, lty=lty, col=col, label.curves=FALSE,
                    ylim=ylim, ...)
  endPlot()
  for(w in c(panel,openPanel)){
    figureCaption = paste("Distribution of time until dropout from",what)
    putFig(w, w, figureCaption,
           if(w==panel) paste(figureCaption,". \\protect\\treatkey", sep="")
           else figureCaption, append=FALSE)
  }
  invisible()
}
