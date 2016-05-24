#' @export
rev.timedom = function(x){
  x$lags = rev(x$lags)
  x
}

#' @export
rev.freqdom = function(x){
  x$freq = rev(x$freq)
  x
}
