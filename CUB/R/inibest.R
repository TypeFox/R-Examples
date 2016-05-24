#' @title Preliminary estimators for CUB models without covariates
#' @description Compute preliminary parameter estimates of a CUB model without covariates for given ordinal
#'  responses. These preliminary estimators are used within the package code to start the E-M algorithm.
#' @aliases inibest
#' @usage inibest(m, freq)
#' @param m Number of ordinal categories
#' @param freq Vector of the absolute frequencies of given ordinal responses
#' @export inibest
#' @return A vector \eqn{(\pi,\xi)} of the initial parameter estimates for a CUB model without covariates,
#'  given the absolute frequency distribution of ordinal responses
#' @seealso \code{\link{inibestgama}}
#' @references Iannario M. (2008). Selecting feeling covariates in rating surveys, \emph{Rivista di Statistica Applicata}, 
#' \bold{20}, 103--116 \cr
#'Iannario M. (2009). A comparison of preliminary estimators in a class of ordinal data models, 
#' \emph{Statistica & Applicazioni}, \bold{VII}, 25--44 \cr
#'  Iannario M. (2012). Preliminary estimators for a mixture model of ordinal data, 
#'  \emph{Advances in Data Analysis and Classification}, \bold{6}, 163--184
#' @keywords htest utilities
#' @examples
#' m<-9
#' freq<-c(10,24,28,36,50,43,23,12,5)
#' estim<-inibest(m, freq) 
#' pai<-estim[1]
#' csi<-estim[2]



inibest <-
function(m,freq){
  freq<-freq/sum(freq)
  csi<-1+(0.5-which.max(freq))/m 
  ppp<-probbit(m,csi)
  pai<-sqrt((sum(freq^2)-1/m)/(sum(ppp^2)-1/m))
  pai<-min(pai,0.99)
  return(c(pai,csi))
}
