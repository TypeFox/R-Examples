#' @title Plot a Discovery Curve
#' @description Plot a species discovery curve.
#' 
#' @param x result of a call to \code{discovery.curve}. 
#' @param col color of confidence interval polygon and line denoting 
#'   \code{s.est}.
#' @param lwd line widths.
#' @param xlab,ylab labels of x and y axes. Only used if \code{add = FALSE}.
#' @param add logical. If TRUE, polygon and lines are added to the current plot.
#' @param ... other arguments passed to plot (ignored).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @references Colwell, R.K., A. Chao, N.J. Gotelli, S.-Y. Lin, C.X. Mao, 
#'   R.L. Chazdon, and J.T. Longino. 2012. Models and estimators linking 
#'   individual-based and sample-based rarefaction, extrapolation and 
#'   comparison of assemblages. Journal of Plant Ecology 5(1):3-21.
#'   
#' @seealso \code{\link{discovery.curve}}
#'
#' @examples
#' data(osa.old.growth)
#' f <- expand.freqs(osa.old.growth)
#' d <- discovery.curve(f, f0.func = Chao1, max.x = 1200)
#' plot(d)
#' 
#' @importFrom graphics plot.new plot.window axis mtext polygon lines points
#' @export plot.discovery.curve
#' @export
#' 
plot.discovery.curve <- function(x, col = "darksalmon", lwd = 2, 
                                 xlab = "# Samples", ylab = "n", 
                                 add = FALSE, ...) {
  s.obs <- x$f.stats["s.obs"]
  s.est <- x$f.stats["s.obs"] + x$f.stats["f0"]
  n <- x$f.stats["n"]
  
  if(!add) {
    plot.new()
    plot.window(xlim = range(pretty(x$s.ind[, "m"])), 
                ylim = range(pretty(c(s.est, x$s.ind.ci))))
    axis(1, pos = 0, lwd = 2)
    axis(2, pos = 0, lwd = 2)
    mtext(xlab, 1, line = 2)
    mtext(ylab, 2, line = 2)
  }
  
  polygon(x$ci.poly, col = col, border = NA)
  lines(x$rarefact.line, lwd = lwd, lty = "solid")
  lines(x$extrap.line, lwd = lwd, lty = "dashed")
  points(n, s.obs, pch = 19)
  
  invisible(NULL)
}