#' Checks if a given object is a time domain filter or process.
#'
#' @title Check if a given object is a time domain object
#' @param X an object
#' @return boolean
#' @export 
is.timedom = function (X){
  inClass(X,'timedom')
}
