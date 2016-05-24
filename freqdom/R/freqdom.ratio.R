#' For given spectral densities \code{S} and \code{SC} computes ratio \code{SC \%*\% S^\{-1\}} at each frequency.
#' 
#' @title Compute a ratio of two spectral densities 
#' @param SC first spectral density
#' @param S second spectral density (square)
#' @param n number of observations used for estimation - the precision of inversion is calculated using this parameter
#' @param Kconst used for heuristing in \code{\link{reg.dim.est}}
#' @param K for inversion as in \code{\link{freqdom.inverse}}
#' @return Frequency Domain Operator object
#' @seealso \code{\link{reg.dim.est}}, \code{\link{freqdom.product}}, \code{\link{freqdom.inverse}}
#' @export
#' @examples
#' n = 100
#' X = rar(n)
#' Y = rar(n)
#' SYX = spectral.density(X,Y)
#' SXX = spectral.density(X)
#' R = freqdom.ratio(SYX,SXX,n)
freqdom.ratio = function(SC,S,n=NULL,K=NULL, Kconst = NULL){
  Sinv = freqdom.inverse(S,n=n,K=K,Kconst=Kconst)
  SC %*% Sinv
}