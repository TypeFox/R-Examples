#' @title Main function for CUSH models 
#' @description Main function to estimate and validate a CUSH model for ordinal responses, with or without covariates
#'  to explain the shelter effect.
#' @aliases CUSH
#' @usage CUSH(ordinal, shelter, m=get('m',envir=.GlobalEnv), X = 0, makeplot = TRUE)
#' @param ordinal Vector of ordinal responses
#' @param shelter Category corresponding to the shelter choice
#' @param m Number of ordinal categories (if omitted, it will be assigned to the number of categories
#'  specified in the global environment)
#' @param X Matrix of selected covariates for explaining the shelter effect. If omitted (default), no covariate 
#' is included in the model
#' @param makeplot Logical: if TRUE (default) and if no covariate is included in the model,
#' the algorithm returns a graphical plot comparing fitted probabilities and observed relative frequencies, 
#' and a plot of the log-likelihood function at the final estimate, compared with the log-likelihood values of 
#' the saturated and the uniform models. If only one explicative dichotomous variable is included in the model, 
#' then the function returns a graphical plot comparing the distributions of the responses conditioned to the 
#' value of the covariate
#' @export CUSH
#' @return An object of the class "CUSH" is a list containing the following results: 
#' \item{estimates}{Maximum likelihood parameters estimates}
#' \item{loglik}{Log-likelihood function at the final estimates}
#' \item{varmat}{Variance-covariance matrix of final estimates (if X=0, it returns the square of the estimated standard error 
#' for the shelter parameter \eqn{\delta})}
#' \item{BIC}{BIC index for the estimated model}
#' @details The estimation procedure is not iterative, so a null result for CUSH$niter is produced.
#' The optimization procedure is run via "optim". If covariates are included, the variance-covariance matrix 
#' is computed as the inverse of the returned numerically differentiated Hessian matrix (option: hessian=TRUE
#'  as argument for "optim"). If not positive definite, it returns a warning message and produces a matrix 
#'  with NA entries.
#' @references 
#' Capecchi S. and Piccolo D. (2015). Dealing with heterogeneity/uncertainty in sample survey with ordinal data, 
#' \emph{IFCS Proceedings, University of Bologna} \cr
#' Capecchi S and Iannario M. (2015). Gini heterogeneity index for detecting uncertainty in ordinal data surveys,
#'  \emph{SIS Proceedings, Treviso - Ca'Foscari, University of Venice}
#' @seealso \code{\link{cushforsim}}, \code{\link{loglikCUSH}}
#' @keywords models
#' @examples
#' data(relgoods)
#' dog<-na.omit(relgoods[,49])
#' m<-10
#' shelter<-1
#' model<-CUSH(dog,shelter=shelter)
#' delta<-model$estimates # ML estimates of delta
#' maxlik<-model$loglik   # Log-likelihood at ML estimates
#' errst<-model$varmat    # Square of estimated standard error of delta
#' BIC<-model$BIC
#' ###############################################
#' ### CUSH model with covariates
#' music<-relgoods[,47]
#' shelter<-1
#' cov<-relgoods[,12]
#' nona<-na.omit(cbind(music,cov))
#' ordinal<-nona[,1]
#' smoking<-nona[,2]
#' modelcov<-CUSH(ordinal,shelter=shelter,X=smoking)
#' omega<-modelcov$estimates
#' maxlik<-modelcov$loglik
#' varmat<-modelcov$varmat
#' BIC<-modelcov$BIC


CUSH <-
function(ordinal,shelter,m=get('m',envir=.GlobalEnv),X=0,makeplot=TRUE){  
  
  rx<-NROW(X)
  if(rx==1){ 
    cush00(m,ordinal,shelter,makeplot)
  }
  else{
    cushcov(m,ordinal,X,shelter,makeplot)
  }
}
