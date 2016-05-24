#Author: Antonio, Fabio Di Narzo. Last Modified 30 March 2011


#'Bivariate time series plots
#'
#'Bivariate time series plots: scatterplots, directed lines and kernel density
#'estimations using functions in the \pkg{sm} package.
#'
#'Bivariate time series plots: scatterplots, directed lines and kernel density
#'and regression functions estimations using functions in the package \pkg{sm}.
#'In particular, for kernel density estimation \code{\link[sm]{sm.density}} is
#'used, with smoothing parameter \code{h} defaulting to
#'\code{\link[sm]{hnorm}}. For kernel regression,
#'\code{\link[sm]{sm.regression}} is used.
#'
#'@param x time series
#'@param lag time lag
#'@param h kernel window (useful only for kernel estimations)
#'@param type type of plot: contour levels, perspective plots, image, directed
#'lines, points or points with superposed kernel regression
#'@return None. Plots are produced on the default graphical device.
#'@author Wrappers to \pkg{sm} by Antonio, Fabio Di Narzo
#'@seealso For finer control on density estimation, consider using directly
#'\code{\link[sm]{sm.density}} and, especially, \code{\link[sm]{sm.ts.pdf}}
#'from package \pkg{sm}.
#'@keywords ts
#'@export
#'@examples
#'
#'x <- log10(lynx)
#'autopairs(x, lag=2, type="lines")
#'
autopairs <- function(x, lag=1, h,
                      type=c("levels","persp","image","lines","points","regression")) {
  panel <- list(levels = function()  sm::sm.density(X, h=rep(h,2), xlab=xlab, ylab=ylab, main="density", display="slice"),
		persp = function() sm::sm.density(X, h=rep(h,2), xlab=xlab, ylab=ylab, main="density", display="persp"),
		image = function() sm::sm.density(X, h=rep(h,2), xlab=xlab, ylab=ylab, main="density", display="image"),
		lines = function() plot(X, xlab=xlab, ylab=ylab, main="lines", type="l"),
		points = function() plot(X, xlab=xlab, ylab=ylab, main="scatter"),
		regression = function() sm::sm.regression(X[,1], X[,2], h=h, xlab=xlab, ylab=ylab, main="regression", ask=FALSE))
  lags <- c(-lag, 0)
  X <- embedd(x, lags=lags)
  xlab <- paste("lag",lag)
  ylab <- paste("lag",0)
  type <- match.arg(type)
  if(missing(h)) {
    h <- sm::hnorm(X)[1]
  }
  panel[[type]]()
}
