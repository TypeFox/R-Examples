#' Creates a frequency domain operator
#'
#' @title Create a frequency domain operator
#' @param X array of evaluations
#' @param freq frequencies at which the object was evaluated
#' @return Frequency domain object
#' @export
freqdom = function (X,freq)
{
  if (!is.array(X) || length(dim(X)) != 3)
    stop("X must be an array of evaluations")
  
  res = list()
  res$operators = X
  res$freq = freq
  class(res) = 'freqdom'
  res
}
