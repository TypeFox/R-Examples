#' Plot criterion-related profile
#'
#' Plots the criterion-related level and pattern profiles for each observation
#'
#' @method plot critpat
#' @importFrom graphics plot abline lines text
#' @export
#' @param x \code{critpat} object resulting from \code{cpa}
#' @param ... additional arguments affecting the plot produced.
#' @seealso \code{\link{cpa}}
plot.critpat <- function(x, ...){
  if(is.null(x$xc) == TRUE)
  {
    plot(x=x$b[[1]],type="b",col = "black", ylim = c(-1,1), ylab = "Estimated Parameters")
    abline(a=0,b=0)
    lines(x=x$b[[2]],type = "b", col = "red")
    text(x$b[[1]],labels(x$b[[1]]), cex = 1, pos = 1)
  }
  if(is.null(x$xc) == FALSE){
  plot(x$xc, type="b",ylab="Criterion Pattern Score")
  abline(a=0,b=0)
  text(x$xc,labels(x$xc), cex=1, pos=3)
  }
}