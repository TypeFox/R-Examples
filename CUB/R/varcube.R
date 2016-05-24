#' @title Variance of CUBE models without covariates
#' @description Compute the variance of a CUBE model without covariates.
#' @aliases varcube
#' @usage varcube(m, pai, csi, phi)
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter 
#' @param phi Overdispersion parameter  
#' @export varcube
#' @seealso  \code{\link{CUBE}}, \code{\link{probcube}}, \code{\link{expcube}}
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in 
#' Ordinal Data,  \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 
#' 771--786
#' @keywords distribution
#' @examples 
#' m<-7
#' pai<-0.8
#' csi<-0.2
#' phi<-0.05
#' varianceCUBE<-varcube(m, pai, csi, phi)

varcube <-function(m,pai,csi,phi){
  odeffect<-pai*csi*(1-csi)*(m-1)*(m-2)*phi/(1+phi)
  variance<- varcub00(m,pai,csi) + odeffect
  return(variance)
}
