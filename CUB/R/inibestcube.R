#' @title Naive estimates for CUBE models without covariates
#' @description Compute \emph{naive} parameter estimates of a CUBE model without covariates for given ordinal responses. 
#' These preliminary estimators are used within the package code to start the E-M algorithm.
#' @aliases inibestcube
#' @usage inibestcube(m,ordinal)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @export inibestcube
#' @return A vector \eqn{(\pi, \xi ,\phi)} of parameter estimates of a CUBE model without covariates
#' @keywords htest utilities
#' @seealso \code{\link{inibestcubecov}}, \code{\link{inibestcubecsi}}
#' @examples
#' data(relgoods)
#' m<-10
#' ordinal<-na.omit(relgoods[,44])
#' estim<-inibestcube(m,ordinal)     # Preliminary estimates: (pai,csi,phi)


inibestcube <-
function(m,ordinal){
  freq<-tabulate(ordinal,nbins=m)  
  inipaicsi<-inibest(m,freq)
  pai<-inipaicsi[1];csi<-inipaicsi[2];
  aver<-mean(ordinal) 
  varcamp<-mean(ordinal^2)-aver^2
  varcub<-varcub00(m,pai,csi) 
  phist<-min(max((varcub-varcamp)/(-pai*csi*(1-csi)*(m-1)*(m-2)-varcub+varcamp),0.01),0.5)
  initial<-c(pai,csi,phist)
  return(initial)
}
