#' Plot the parameters and its confidence intervals
#' \code{summary} Plot the parameters and its confidence intervals
#' @param object An object for which a summary is desired.
#' @param ... Additional arguments affecting the summary produced.
#' @export
summary.quickpsy <- function(object,...)
{
  print(object$par)
  print(object$parci)
}
