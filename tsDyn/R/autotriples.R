#Author: Antonio, Fabio Di Narzo. Last Modified 30 March 2011


#'Trivariate time series plots
#'
#'Trivariate time series plots: kernel autoregression using functions in the
#'\pkg{sm} package
#'
#'This function displays trivariate time series plots, i.e. kernel regression
#'of \eqn{x[t-lags[1]], x[t-lags[2]]}{x_{t-l_1}, x_{t-l_2}} against
#'\eqn{x[t]}{x_t} using functions in the package \pkg{sm}.  In particular,
#'\code{\link[sm]{sm.regression}} is used, with smoothing parameter defaulting
#'to \code{\link[sm]{hnorm}(x)}.
#'
#'@param x time series
#'@param lags vector of regressors lags
#'@param h kernel window
#'@param type type of plot: contour levels, perspective plots, image
#'@return None. Plots are produced on the default graphical device.
#'@author Wrappers to \pkg{sm} by Antonio, Fabio Di Narzo
#'@seealso For finer control on kernel regression, consider using directly
#'\code{\link[sm]{sm.regression}} and, especially,
#'\code{\link[sm]{sm.autoregression}} in package \code{\link[sm]{sm}}.
#'@keywords ts
#'@export
#'@examples
#'
#'autotriples(log(lynx))
#'autotriples(log(lynx), type="persp")
#'autotriples(log(lynx), type="image")
#'
autotriples <- function(x, lags=1:2, h, type=c("levels","persp","image", "lines", "points")) {
  
  panel <- list(levels = function(x) contour(x, xlab=xlab, ylab=ylab),
		persp = function(x) persp(x, xlab=xlab, ylab=ylab, zlab=zlab),
		image = function(x) image(x, xlab=xlab, ylab=ylab),
		lines = function(x) scatterplot3d::scatterplot3d(X, xlab=xlab, ylab=ylab, zlab=zlab, main="directed lines", type="l"),
		points = function(x) scatterplot3d::scatterplot3d(X, xlab=xlab, ylab=ylab, zlab=zlab, main="cloud", pch=1))
  type <- match.arg(type)
  X <- embedd(x, lags=c(-lags,0))
  if(missing(h)) 
    h <- sm::hnorm(X[,1])
  xlab <- paste("lag",lags[1])
  ylab <- paste("lag",lags[2])
  zlab <- "lag 0"
  mod <- sm::sm.regression(X[,1:2], X[,3], h=rep(h,2), display="none")
  panel[[type]](mod$estimate)
  invisible(NULL)
}
