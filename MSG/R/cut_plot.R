#' Cut the points in a scatter plot into groups according to x-axis
#'
#' This function can categorize the variable on the x-axis into groups and plot
#' the mean values of y. The purpose is to show the arbitrariness of the
#' discretization of data.
#' @param x the x variable
#' @param y the y variable
#' @param breaks the breaks to cut the x variable
#' @param ... other arguments to be passed to
#'   \code{\link[graphics]{plot.default}}
#' @param pch.cut the point symbol to denote the mean values of y
#' @return NULL
#' @export
#' @author Yihui Xie <\url{http://yihui.name}>
#' @examples x=rnorm(100); y=rnorm(100)
#' cut_plot(x,y,seq(min(x),max(x),length=5))
cut_plot = function(x, y, breaks, ..., pch.cut = 20) {
  plot(x, y, ...)
  x1 = cut(x, breaks)
  mid = (breaks[-length(breaks)] + breaks[-1])/2
  abline(v = breaks, lty = 2)
  points(mid,tapply(y, x1, mean), type = 'o', pch = pch.cut)
}
