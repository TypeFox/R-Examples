#' @title Grid-based preliminary parameter estimates for CUB models
#' @description Compute the log-likelihood function of a CUB model with parameter vector \eqn{(\pi, \xi)} ranging in
#' the Cartesian product between \eqn{x} and \eqn{y}, for a given absolute frequency distribution.
#' @aliases inigrid
#' @usage inigrid(m, freq, x, y)
#' @param m Number of ordinal categories
#' @param freq Vector of length \eqn{m} of the absolute frequency distribution
#' @param x A set of values to assign to the uncertainty parameter \eqn{\pi}
#' @param y A set of values to assign to the feeling parameter \eqn{\xi}
#' @export inigrid
#' @return It returns the parameter vector corresponding to the maximum value of the log-likelihood 
#' for a CUB model fitting ordinal responses with given frequencies.
#' @seealso \code{\link{inibest}} 
#' @keywords htest utilities
#' @examples
#' m<-9
#' x<-c(0.1,0.4,0.6,0.8)
#' y<-c(0.2, 0.5,0.7)
#' freq<-c(10,24,28,36,50,43,23,12,5)
#' ini<-inigrid(m,freq,x,y)
#' pai<-ini[1]
#' csi<-ini[2]


inigrid <-
function(m,freq,x,y){
  listap<-expand.grid(x,y)
  quanti<-NROW(listap)
  loglik<-rep(NA,quanti)
  for(j in 1:quanti){
    pai<-listap[j,1]; csi<-listap[j,2];
    loglik[j]<-loglikcub00(m,freq,pai,csi)
  }
  indice<-which.max(loglik)
  return(listap[indice,])
}
