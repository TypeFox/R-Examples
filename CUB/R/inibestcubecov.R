#' @title Preliminary parameter estimates for CUBE models with covariates
#' @description Compute preliminary parameter estimates for a CUBE model with covariates for all the three parameters. 
#' These estimates are set as initial values to start the E-M algorithm within maximum likelihood estimation.
#' @aliases inibestcubecov
#' @usage inibestcubecov(m, ordinal, Y, W, Z)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of selected covariates to explain the uncertainty parameter
#' @param W Matrix of selected covariates to explain the feeling parameter
#' @param Z Matrix of selected covariates to explain the overdispersion parameter
#' @export inibestcubecov
#' @return A vector (inibet, inigama, inialpha) of preliminary estimates of parameter vectors for 
#'  \eqn{\pi = \pi(\beta)}, \eqn{\xi=\xi(\gamma)}, \eqn{\phi=\phi(\alpha)}, respectively, of a CUBE model with covariates for all the three
#'   parameters. In details, inibet, inigama and inialpha have length equal to NCOL(Y)+1, NCOL(W)+1 and
#'   NCOL(Z)+1, respectively, to account for an intercept term for each component 
#' @keywords htest utilities
#' @seealso \code{\link{inibestcube}}, \code{\link{inibestcubecsi}}, \code{\link{inibestgama}}
#' @examples
#' data(relgoods)
#' m<-10
#' covpai<-relgoods[,2]
#' covcsi<-relgoods[,6]
#' covphi<-relgoods[,7]
#' ordinal<-relgoods[,37]
#' nona<-na.omit(cbind(ordinal,covpai,covcsi,covphi))  # Omitting missing values
#' ordinal<-nona[,1]
#' Y<-nona[,2]
#' W<-nona[,3]
#' Z<-nona[,4]
#' ini<-inibestcubecov(m, ordinal, Y, W, Z)
#' p<-NCOL(Y)
#' q<-NCOL(W)
#' inibet<-ini[1:(p+1)]               # Preliminary estimates for uncertainty 
#' inigama<-ini[(p+2):(p+q+2)]        # Preliminary estimates for feeling 
#' inialpha<-ini[(p+q+3):length(ini)] # Preliminary estimates for overdispersion


inibestcubecov <-
function(m,ordinal,Y,W,Z){
  q<-NCOL(Y)
  p<-NCOL(W)
  v<-NCOL(Z)
  inigama<-inibestgama(m,ordinal,W)
  inicube<-inibestcube(m,ordinal)
  pai<-inicube[1]
  bet0<-log(pai/(1-pai))
  inibet<-c(bet0,rep(0.1,q))
  alpha0<-log(0.1)
  inialpha<-c(alpha0,rep(0.1,v))
  return(c(inibet,inigama,inialpha))
}
