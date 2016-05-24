# get one lag: A_lag = \int R e^{itlag}
inv.fourier.one = function(R,lag){
  A = array(0,dim(R$operators)[2:3])
  A[,] = 0
  for (theta in 1:length(R$freq))
    A = A + R$operators[theta,,] * exp(R$freq[theta]*1i*lag)
  A / length(R$freq)
}

#' Inverse Fourier transform of a given Frequency Domain Operator.
#'
#' @title Inverse Fourier transform for operators
#' @param R the operator to invert
#' @param lags lags to compute
#' @return Time domain operator serie object
#' @export
#' @seealso \code{\link{fourier.transform}}
#' @examples
#' n = 100
#' X = rar(n)
#' Y = rar(n)
#' #estimate regressors in model $Y_t = \sum_{i\in Z} A_i X_{t-i}$
#' SYX = spectral.density(Y, X)
#' SXX = spectral.density(X)
#' R = freqdom.ratio(SYX,SXX, n)
#' A = invfourier(R) 
invfourier = function(R,lags=0:0){
  if (!is.freqdom(R))
    stop("R must be a freqdom object")
  if (!is.numeric(lags) || !is.vector(lags))
    stop("lags must be a vector of integers")
  
  H = length(lags)
  A = array(0, dim=c(H, dim(R$operators)[2:3]))
  for (h in 1:H)
    A[h,,] = inv.fourier.one(R,lags[h])
  
  timedom(Re(A),lags)
}
