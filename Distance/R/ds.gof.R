#' Goodness of fit tests for distance sampling models
#'
#' Chi-square, Kolmogorov-Smirnov and Cramer-von Mises goodness of fit tests for detection function models.
#'
#' @export
#' @param model fitted model object
#' @param breaks Cutpoints to use for binning data
#' @param nc Number of distance classes
#' @param qq Flag to indicate whether quantile-quantile plot is desired
#' @param \dots Graphics parameters to pass into qqplot function
#' @return List of test results and a plot.
#' @author David L Miller
#' @seealso qqplot.ddf ddf.gof
#' @keywords utility
ds.gof <- function(model, breaks=NULL, nc=NULL, qq=TRUE,...){
  return(ddf.gof(model$ddf, breaks=breaks, nc=nc, qq=qq,...))
}
