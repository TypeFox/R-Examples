#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2006-03-04 22:04:02 +0100 (sab, 04 mar 2006) $


#'Interactive trivariate time series plots
#'
#'Interactive trivariate time series plots
#'
#'This function displays interactive trivariate time series plots
#'\code{x[t-lags[1]], x[t-lags[2]]} against \code{x[t]} using the interactive
#'\code{\link[rgl]{rgl}} device.
#'
#'@param x time series
#'@param lags vector of regressors lags
#'@param type type of plot: contour levels, perspective plots, image
#'@return None. A plot is produced on the current \code{rgl} device.
#'@author Wrapper to 'sm' and GUI by Antonio, Fabio Di Narzo
#'@seealso \code{\link{autotriples}} for 3d visualization via
#'\pkg{scatterplot3d} package and for kernel post-processing of the cloud for
#'nonparametric autoregression functions estimates.
#'@keywords ts
#'@export
#'@examples
#'
#'if(interactive())
#' autotriples.rgl(log(lynx))
#'
autotriples.rgl <- function(x, lags=1:2, type=c("lines","points")) {
	type <- match.arg(type)
	X <- embedd(x, lags=c(-lags,0))
	rgl::rgl.clear()
	if(type=="lines")
	  rgl::rgl.linestrips(X[,1],X[,2],X[,3])
	else if (type=="points")
	  rgl::rgl.points(X[,1],X[,2],X[,3])
        invisible(NULL)
}
