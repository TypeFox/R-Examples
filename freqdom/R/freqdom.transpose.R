#' Transpose the object at each frequency or lag
#'  
#' @title Transpose pointwise timedom or freqdom object
#' @param x freqdom or timedom object
#' @return object of the same type as x
#' @export
freqdom.transpose = function(x){
  lags = freqdom.lags(x)
  for (i in 1:length(lags))
    x$operators[i,,] = t(x$operators[i,,])
  x
}

#' @S3method t freqdom
t.freqdom = freqdom.transpose

#' @S3method t timedom
t.timedom = freqdom.transpose
