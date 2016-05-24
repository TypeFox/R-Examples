#' Compute Fourier transform of a given Time Domain Operator
#'
#' @title Compute fourier transform of given series of operators 
#' @param A time domain operator series
#' @param freq frequencies on which the transfom should be computed
#' @seealso \code{\link{invfourier}}
#' @examples
#' X = rar(100,d=2)
#' C = cov.structure(X)
#' F = fourier.transform(C) # a simple spectral density estimator
#' Cinv = invfourier(F)
#' @export
fourier.transform = function(A,freq=NULL){
  if (is.vector(A) || is.matrix(A))
    A = timedom(A)
  if (!is.timedom(A))
    stop("A must be a time domain object")
  
  if (is.null(freq))
    freq = pi*-100:100/100
  thetas = freq

  if (!is.vector(thetas) || !is.numeric(thetas) ||
        max(thetas) > pi + 0.000001 || min(thetas) < - pi - 0.000001)
    stop("freq must be a vector of real number from [-pi,pi] intercval")
  
  D = dim(A$operators)
  nbasisX = D[2]
  nbasisY = D[3]
  lags = A$lags
  
  S = array(0,c(length(thetas),nbasisX,nbasisY))
  
  # Compute sum at each frequency
  for (theta in 1:length(thetas))
  {
    # Compute one summand
    H = length(A$lags)
    
    for (h in 1:H){
      lag = lags[h]
      R = exp(-1i*lag*thetas[theta]) * A$operators[h,,]
      S[theta,,] = S[theta,,] + R
    }
    S[theta,,] = (S[theta,,]) #/(max(thetas) - min(thetas)) # TODO: CHECK
  }
  freqdom(S,thetas)
}
