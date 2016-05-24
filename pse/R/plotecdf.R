#' Uncertainty and Sensitivity Plots.
#'
#' The functions listed here are used in uncertainty and sensitivity estimation.
#' 
#' The function \code{plotscatter} produces a series of scatterplots from data.
#' 
#' The function \code{plotecdf} plots the empirical cumulative density function
#'	 from an LHS object	or PLUE object. 
#'
#' The function \code{plotprcc} plots the partial rank correlation coefficient
#'	 from an LHS object	or PLUE object. 
#' 
#'  Finally, the \code{plotcv} function plots the empirical cummulative density function 
#'	(ecdf) of the coefficient of variation of the LHS resulting from a stochastic 
#'	simulation, along with a dotted line representing the coefficient of variation
#'	of the whole result set. See the 'multiple' vignette for examples and interpretation.
#' 
#' The function plotscatter accepts an alternative invocation of \code{plotscatter(obj, res)}
#' in which obj is a data.frame consisting on the data to be plotted on the x axis, and 
#' res is a data.frame consisting on the model results to be plotted on the y axis.
#' 
#' @param obj The LHS or PLUE object containing the simulation results to be plotted.
#'
#'	  NOTICE: plotecdf and plotcv only accept LHS objects! For plotting the likelihood profile
#'	  from a PLUE object, simply use \code{plot(obj)}
#' @param res A data.frame consisting of the model results to be plotted on the y axis, if
#' 'obj' is passed as a data.frame. If 'obj' is an LHS/PLUE object, this parameter is ignored.
#' @param stack If the results is a data.frame with several variables, \code{stack=FALSE} generates
#'	  a series of plots, and \code{stack=TRUE} generates a single plot with the ECDF from
#'	  all variables identified by different colors.
#' @param index.res An optional vector indicating which columns from the results are to be plotted.
#' @param index.data An optional vector with the indices of the data columns to be plotted. 
#' @param col An optional vector indicating the colors to be used.
#' @param xlab,ylab Labels for the x axis (ecdf) or y axis(prcc). 
#'	  The functions use the name provided in the res.names argument from the LHS function if left blank.
#' @param add.lm Boolean. Whether to include a simple linear model on the plots. Defaults to TRUE.
#' @param \dots Additional parameters to be passed to the lower level plotting function.
#' @examples
#' myLHS <- LHS(model=function(x) x[,1]+x[,2]*x[,3], factors=3, N=20, res.names="My Output")
#' plotecdf(myLHS, main="ECDF plot")
#' plotprcc(myLHS, main="PRCC plot")
#' plotscatter(myLHS)
#' @export
#' @rdname plots
#' @import graphics Hmisc
plotecdf <- function (obj, stack=FALSE, index.res =1:get.noutputs(obj), col=index.res, xlab = NULL, ...) {
	if (is.null (xlab)) xlab = obj$res.names
	if (stack) {
		if (length(xlab) > 1) xlab = "obj results"
		dat <- vec(get.results(obj)[,index.res])
		g <- rep(index.res, each=dim(obj$res)[1])
		Ecdf(dat, group=g, col=col, xlab=xlab, ...)
	} else Ecdf(get.results(obj)[,index.res], xlab=xlab, ...)
}

