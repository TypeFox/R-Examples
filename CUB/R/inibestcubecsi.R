#' @title Preliminary estimates of parameters for CUBE models with covariates only for feeling
#' @description Compute preliminary parameter estimates of a CUBE model with covariates only for feeling, given
#'  ordinal responses. These estimates are set as initial values to start the E-M algorithm for such model.
#' @aliases inibestcubecsi
#' @usage inibestcubecsi(m, ordinal, W, starting, maxiter, toler)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param W Matrix of selected covariates to explain the feeling parameter
#' @param starting Starting values for \code{\link{cubeforsim}} simulation routine
#' @param maxiter Maximum number of iterations allowed for running the \code{\link{cubeforsim}} simulation routine 
#' @param toler Fixed error tolerance for final estimates for \code{\link{cubeforsim}} simulation routine
#' @export inibestcubecsi
#' @details It invokes \code{\link{cubeforsim}} to obtain preliminary estimates
#'  for the uncertainty and the overdispersion parameters. As to the feeling component, it considers the
#'   nested CUB model with covariates and calls \code{\link{inibestgama}} to derive initial estimates for the coefficients
#'   of the selected covariates.
#' @return A vector (pai, gamaest, phi), where pai is the initial estimate for the uncertainty parameter, gamaest
#'  is the vector of initial estimates for the feeling component (including an intercept term in the first entry),
#'   and phi is the initial estimate for the overdispersion parameter
#' @keywords htest utilities
#' @seealso \code{\link{inibestcube}}, \code{\link{inibestcubecov}}, \code{\link{inibestgama}}, \code{\link{cubeforsim}}
#' @examples
#' data(relgoods)
#' isnacov<-which(is.na(relgoods[,2]))
#' isnaord<-which(is.na(relgoods[,37]))
#' unina<-union(isnacov,isnaord);
#' newdati<-relgoods[-unina,]
#' ordinal<-newdati[,37]
#' W<-newdati[,2]
#' m<-10
#' starting<-rep(0.1, 3)
#' ini<-inibestcubecsi(m, ordinal, W, starting, maxiter=100, toler=1e-3)
#' nparam<-length(ini)
#' pai<-ini[1]                 # Preliminary estimates for uncertainty component
#' gamaest<-ini[2:(nparam-1)]  # Preliminary estimates for coefficients of feeling covariates
#' phi<-ini[nparam]            # Preliminary estimates for overdispersion component


inibestcubecsi <-
function(m,ordinal,W,starting,maxiter,toler){
  gamaest<-inibestgama(m,ordinal,W)
  stimacube<-cubeforsim(m,ordinal,starting,maxiter,toler,expinform=FALSE)
  param<-stimacube$estimates 
  elle<-length(param)
  pai<-param[1]; phi<-param[elle];
  iniest<-c(pai,gamaest,phi)
  return(iniest)
}
