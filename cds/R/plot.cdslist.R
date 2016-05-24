#' Plot a \code{cdslist} Object
#' 
#' Create a scree plot and bubble plots for all elements in a \code{cdslist} object.
#' 
#' @param x An object of class \code{cdslist}.
#' @param which The which argument passed to \code{\link{plot.cds}}.
#' @param \dots Additional arguments passed to \code{\link{plot.cds}}.
#' @method plot cdslist
#' @export
plot.cdslist <- function(x, which = 2L, ...){
  loss <- sapply(x, "[[", "minloss")
  K <- sapply(x, "[[", "K")
  plot(K, loss, type = "b", pch = 16, main = "Scree Plot", xlab = "K", ylab = "Loss")
  invisible(lapply(x, plot, which = which, ...))
}