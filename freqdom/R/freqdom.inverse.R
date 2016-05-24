#' For given Frequency Domain Operator \code{S} computes ratio \code{S^\{-1\}} at each frequency.
#' 
#' @title Compute an inverse of a given Frequency Domain Operator
#' @param S first spectral density
#' @param n number of observations used for estimation - the precision of inversion is calculated using this parameter
#' @param K how many directions should be inverted (as in \code{\link{pseudoinverse}})
#' @param Kconst constant for heuristic (as in \code{\link{reg.dim.est}})
#' @return Frequency Domain Operator object
#' @seealso \code{\link{freqdom.product}}, \code{\link{freqdom.ratio}}
#' @export
freqdom.inverse = function(S,n=NULL,Kconst = NULL, K=NULL){
  if (!is.freqdom(S))
    stop("S must be a freqdom object")
  if (!is.null(K) && !is.positiveint(K+1) && is.vector(K))
    stop("K must be a nonnegative integer")
  
  # Compute number of directions to estimate
  E = freqdom.eigen(S)
  
  if (is.null(K)){
    th = 1/n
    
    K = sum(freqdom.inflambdas(E) / freqdom.inflambdas(E)[1] >= th)
    debug.trace(paste("n =",n,", we have chosen K",K))
  }
  if (K <= 1 && !is.vector(K))
    K = 1

  R = S
  for (theta in 1:length(S$freq)){
    th = 0
    # Instead of inversing S we just take several eigendirections
    
    D = dim(S$operators)
    A = matrix(S$operators[theta,,],D[2],D[3])
    
    curK = K
    if (is.vector(K) && length(K) > 1)
      curK = K[theta]
    R$operators[theta,,] = pseudoinverse(A,curK)
  }
  R
}
