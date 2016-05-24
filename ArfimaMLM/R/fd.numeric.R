#' @export
#' @importFrom fracdiff diffseries
fd.numeric <-
function(x, dval, ...) {
  if(dval!=0) {
    if (dval>1) dval <- 1
    x <- diffseries(x, dval)
  }
  out <- list(series=x,estimator="external",value=dval)
  out
}
