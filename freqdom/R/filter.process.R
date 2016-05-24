#' Compute
#' \deqn{Y_t = \sum_{k=-q}^p A_k X_{t-k}}
#' where \eqn{X_t} is a stationary multivariate time series and \eqn{(A_k)_{-q \leq k \leq p}} is a filter.
#' Note that this procedure can also be used for deconvolution 
#' \deqn{\tilde X_t = \sum_{k=-q}^p X_{t+k}} A_k'
#' but therefore in place of \code{A} one should it's transposition and inversed time, i.e.
#' \code{t(rev(A))}.
#'
#' @title Generate a linear process from 
#' @param X process
#' @param A time-domain operator series
#' @return Multivariate linear process
#' @export
filter.process = function(X, A){
  if (!is.timedom(A))
    stop("A must be timedom object")
  if (!is.matrix(X))
    stop("X must be a multivariate time series (a matrix of observations)")
  
  n = dim(X)[1]
  nbasis = dim(A$operators)[2]
  Y = matrix(0,n,nbasis)
  
  for (component in 1:nbasis)
  {
    IP = A$operators[,component,] %*% t(X)
    for (i in 1:length(A$lags)){
      lag = A$lags[i]
      IP[i,] = colshift(IP[i,],lag)
    }
    Y[,component] = colSums(IP)
    
  }
  Y
}

colshift = function(col,lag){
  n = length(col)
  if (abs(lag) >= n)
    rep(0,n)
  else if (lag<0)
    c(col[-(1:(-lag))],rep(0,(-lag)))
  else if (lag>0)
    c(rep(0,lag),col[(lag+1):n - lag])
  else
    col
}


#' Convolution of a process X with an operator A. As in \code{\link{filter.process}}
#'  
#' @title Convolution of a process X with an operator A.
#' @param A \code{timedom} operators
#' @param X \code{timedom} multivariate time series
#' @return Convoluted series of the same type as X
#' @export
`%c%` <- function(A, X){
  filter.process(X, A)
}
