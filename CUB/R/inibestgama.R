#' @title Preliminary parameter estimates of a CUB model with covariates for feeling
#' @description Compute preliminary parameter estimates for CUB models 
#' with covariates for feeling, given ordinal responses.
#' These estimates are set as initial values for parameters to start the E-M algorithm.
#' @aliases inibestgama
#' @usage inibestgama(m, ordinal, W)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param W Matrix of selected covariates for explaining the feeling component
#' @export inibestgama
#' @return A vector of length equal to NCOL(W)+1, whose entries are the preliminary estimates
#'  of the parameters for the feeling component, including an intercept term as first entry
#' @references Iannario M. (2008). Selecting feeling covariates in rating surveys, 
#' \emph{Rivista di Statistica Applicata}, \bold{20}, 103--116 \cr
#' Iannario M. (2009). A comparison of preliminary estimators in a class of ordinal data models,
#' \emph{Statistica and Applicazioni}, \bold{VII}, 25--44 \cr
#' Iannario M. (2012).  Preliminary estimators for a mixture model of ordinal data, 
#' \emph{Advances in Data Analysis and Classification}, \bold{6}, 163--184
#' @seealso \code{\link{inibest}}, \code{\link{inibestcubecsi}}
#' @keywords htest utilities
#' @examples
#' data(univer)
#' m<-7
#' ordinal<-univer[,12]
#' diploma<-univer[,5]
#' ini<-inibestgama(m,ordinal,W=diploma)


inibestgama <-
function(m,ordinal,W){
  WW<-cbind(1,W)                           
  ni<-log((m-ordinal+0.5)/(ordinal-0.5))
  gama<-(solve(t(WW)%*%WW))%*%(t(WW)%*%ni) 
  return(gama)
}
