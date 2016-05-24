inClass = function(X,cls){
  !is.null(oldClass(X)) && oldClass(X) == cls
}

#' Checks if a given object is a frequency domain matrix.
#'
#' @title Check if a given object is a frequency domain matrix
#' @param X an object
#' @return boolean
#' @export 
is.freqdom = function (X){
  inClass(X,'freqdom')
}

is.positiveint = function (n){
  is.numeric(n) && n > 0
}
