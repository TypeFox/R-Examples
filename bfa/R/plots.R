#' @include class.R
NULL

#' Display a biplot
#' @param x A bfa object
#' @param factors Integer vector giving indices of the factors to plot
#' @param ... Additional arguments to biplot; see \code{?biplot}
#' @method biplot bfa
#' @return Shows a biplot
#' @export
biplot.bfa <- function(x, factors=c(1,2), ...) {
  call_args = list(...)
  if(is.null(call_args$xlabs)) call_args$xlabs = x$obslabel
  if(is.null(call_args$ylabs)) call_args$ylabs = x$varlabel
  call_args$x = t(x$post.scores.mean)[,factors]
  call_args$y = x$post.loadings.mean[,factors]
  do.call(biplot, call_args)
}