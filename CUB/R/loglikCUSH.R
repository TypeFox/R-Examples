#' @title Log-likelihood function for CUSH models
#' @aliases loglikCUSH
#' @description  Compute the log-likelihood function for CUSH models with or without covariates 
#' to explain the shelter effect.
#' @usage loglikCUSH(ordinal,m,param,shelter,X=0)
#' @export loglikCUSH
#' @param ordinal Vector of ordinal responses
#' @param m Number of ordinal categories
#' @param param Vector of parameters for the specified CUSH model
#' @param shelter Category corresponding to the shelter choice
#' @param X Matrix of selected covariates to explain the shelter effect (default: no covariate 
#' is included in the model)
#' @details If no covariate is included in the model, then "param" is the estimate of the shelter 
#' parameter (delta), otherwise "param" has length equal to NCOL(X) + 1 to account for an intercept  
#' term (first entry)
#' @seealso  \code{\link{CUSH}}
#' @keywords htest
#' @examples
#' ## Log-likelihood of CUSH model without covariates
#' n<-300
#' m<-7
#' shelter<-2
#' delta<-0.4
#' ordinal<-simcush(n,m,delta,shelter)
#' loglik<-loglikCUSH(ordinal,m,param=delta,shelter)
#' #####################
#' ## Log-likelihood of CUSH model with covariates
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods[,44]
#' cov<-relgoods[,2]
#' nona<-na.omit(cbind(ordinal,cov))
#' ordinal<-nona[,1]
#' cov<-nona[,2]
#' shelter<-1
#' omega<-c(-2.29, 0.62)
#' loglik<-loglikCUSH(ordinal,m,param=omega,shelter,X=cov)
#' 

loglikCUSH<-function(ordinal,m,param,shelter,X=0){
  
  nx<-NROW(X)
  if (nx==1){
    delta<-param
    loglik<-loglikcush00(m,ordinal,delta,shelter)
    
  } else {
    omega<-param
    loglik<-loglikcushcov(m,ordinal,X,omega,shelter)  
  }
  
  return(loglik)
  
}

