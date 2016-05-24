#' @export
print.quickpsy <- function(x,...)
{
  print(x$par)
  print(x$parci)
  if ('thresholds' %in% names(x)) print(x$thresholds)
  if ('thresholdsci' %in% names(x)) print(x$thresholdsci)
}
