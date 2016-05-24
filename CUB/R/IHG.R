#' @title Main function for IHG models 
#' @description Main function to estimate and validate an Inverse Hypergeometric model, without or 
#' with covariates for explaining the preference parameter.
#' @aliases IHG
#' @usage IHG(ordinal, m=get('m',envir=.GlobalEnv), U = 0, makeplot = TRUE)
#' @param ordinal Vector of ordinal responses
#' @param m Number of ordinal categories (if omitted, it will be assigned to the number of categories
#'  specified in the global environment)
#' @param U Matrix of covariates for explaining the preference parameter. If omitted (default), 
#' no covariate is included in the model
#' @param makeplot Logical: if TRUE (default), the algorithm returns a graphical plot comparing fitted probabilities
#'  and observed relative frequencies  for IHG models without covariates. If only one explicative dichotomous 
#'  variable is included in the model, then the function returns a graphical plot comparing the distributions
#'   of the responses conditioned to the value of the covariate
#' @export IHG
#' @return An object of the class "IHG" is a list containing the following results: 
#' \item{estimates}{Maximum likelihood parameters estimates}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates. If no covariate is included in the model, 
#' it returns the square of the estimated standard error for the preference parameter \eqn{\theta}}
#' \item{BIC}{BIC index for the estimated model}
#' @details This is the main function for IHG models (that are nested into CUBE models, see the references below),
#'  calling for the corresponding function whenever covariates are specified. \cr
#' The parameter \eqn{\theta} represents the probability of observing a rating corresponding to the first 
#' category and is therefore a direct measure of preference, attraction, pleasantness toward the investigated item.
#'  This is reason why \eqn{\theta} is customarily referred to as the preference parameter of the IHG model.\cr
#' The estimation procedure is not iterative, so a null result for IHG$niter is produced. \cr
#' The optimization procedure is run via "optim". The variance-covariance matrix (or the estimated standard error for
#'  theta if no covariate is included) is computed as the inverse of the returned numerically differentiated 
#'  Hessian matrix (option: hessian=TRUE as argument for optim). If not positive definite,
#'  it returns a warning message and produces a matrix with NA entries.
#' @references 
#' D'Elia A. (2003). Modelling ranks using the inverse hypergeometric distribution, 
#' \emph{Statistical Modelling: an International Journal}, \bold{3}, 65--78 \cr
#' Iannario M. (2012). CUBE models for interpreting ordered categorical data with overdispersion,
#'  \emph{Quaderni di Statistica}, \bold{14}, 137--140
#' @seealso \code{\link{probihg}},  \code{\link{iniihg}}, \code{\link{loglikIHG}} 
#' @keywords models
#' @examples 
#' \donttest{
#' data(relgoods)
#' m<-10
#' ordinal<-na.omit(relgoods[,41]) 
#' model<-IHG(ordinal)
#' theta<-model$estimates      # ML estimates for the preference parameter theta
#' maxlik<-model$loglik
#' sterr<-model$varmat         # square of estimated standard error for theta
#' BICIHG<-model$BIC
#' #################################
#' ordinal<-relgoods[,41]
#' gender<-relgoods[,9]
#' data<-na.omit(cbind(ordinal,gender))
#' modelcov<-IHG(data[,1],U=data[,2])
#' omega<-modelcov$estimates     #  ML estimates (intercept term: omega[1])
#' maxlik<-modelcov$loglik
#' varmat<-modelcov$varmat
#' BICcov<-modelcov$BIC
#' }

IHG <-
function(ordinal,m=get('m',envir=.GlobalEnv),U=0,makeplot=TRUE){
  
  ru<-NROW(U)
  if (ru==1){
    ihg00(m,ordinal,makeplot)
  } else {
    ihgcov(m,ordinal,U,makeplot)
  }
}
