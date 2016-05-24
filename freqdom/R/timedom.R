#' Creates a time domain operator
#'
#' @title Creates a time domain operator
#' @param X array of operator or matrix of observations
#' @param lags lags on which the operator is defined
#' @return Time domain object
#' @export
timedom = function (X,lags=NULL){
  res = list()

  if (is.vector(X)){
    if (is.null(lags))
      lags = 1:length(X)
    res$operators = array(0,c(length(X),1,1))
    res$operators[,1,1] = X
  }
  else if (is.matrix(X)){
    if (is.null(lags))
      lags = 1:dim(X)[1]

    res$operators = array(0,c(dim(X)[1],dim(X)[2],1))
    res$operators[,,1] = X
  }
  else if (is.array(X) && length(dim(X)) == 3){
    if (is.null(lags))
      lags = 1:dim(X)[3] - 1
    res$operators = X
  }
  else{
    stop("X must be a matrix or an array")
  }
  
  res$lags = lags
  class(res) = 'timedom'
  res
}
