#' Plot the time-dependent ROC curve
#'
#' Plot the ROC curve estimated by \code{tdROC()}.
#' @param x the object returned by \code{tdROC()}.
#' @param lwd user-specified line width. Default is \code{2}.
#' @param xlab user-specified label for x-axis. Default is "\code{1-specificity}".
#' @param ylab user-specified label for y-axis. Default is "\code{sensitivity}".
#' @param xlim user-specified limit for x axis. Default is \code{c(0,1)}.
#' @param ylim user-specified limit for y axis. Default is \code{c(0,1)}.
#' @param col user-specified color for ROC curve. Defualt is "\code{black}".
#' @param main user-specified title for the plot. Default is "\code{ROC curve}"
#' @param abline user-specified reference diagnol line. Default is \code{True}.
#' @param \dots for future methods
#' @return Returns a plot of ROC curve.
#' @importFrom survival survfit
#' @examples
#' library( survival )
#' data( mayo );
#' dat <- mayo[ , c( "time","censor","mayoscore5" )] ;
#' fm <- tdROC( X = dat$mayoscore5, Y = dat$time, delta = dat$censor,
#'              tau = 365*6, span = 0.1, nboot = 0, alpha = 0.05, n.grid = 1000, cut.off = 5:9 )
#' # plot the object "fm" from tdR0C()
#' plot.tdROC( fm ) ;
#'
#' @export

plot.tdROC <- function( x, lwd=2, xlab="1-specificity", ylab="sensitivity", xlim=c(0,1), ylim=c(0,1), col="black", main="ROC curve", abline=T, ...) {
  # Plot the ROC curve estimated by tdROC()
  # Arguments:
  #  -- x: the object returned by tdROC()
  # Return:
  #  -- a plot of ROC curve
  spec <- 1 - x$ROC$spec ;
  sens <- x$ROC$sens ;
  tmp <- order(spec, sens) ;
  x2 <- spec[tmp] ;
  y2 <- sens[tmp] ;
  #windows() ;
  plot( x = x2, y = y2, xlab=xlab, ylab=ylab,
        type="s", lwd = lwd, xlim=xlim, ylim=ylim,
        main=main,
        sub=paste("AUC = ", round(x$AUC$value, 3), " (",
                   round(x$AUC$lower, 3), ", ",
                   round(x$AUC$upper, 3), ")", sep=""),
        col=col
  ) ;
  if (abline==T){
    abline(0, 1, col="gray", lwd=2, lty=2) ;
  }
  invisible(0) ;
}
