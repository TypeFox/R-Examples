#' @name Eriksson
#' @aliases Eriksson
#' @title Eriksson model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Eriksson model. 
#' The function can also return the transformed variable. 
#' The Eriksson model was proposed by Eriksson in 1973.
#' 
#' @usage Eriksson(z,p,mu,sigma,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p probability of direct response
#' @param mu mean of the scramble variable \eqn{S}
#' @param sigma standard deviation of the scramble variable \eqn{S}
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' The randomized response given by the person labelled \eqn{i} is \eqn{y_i} with probability \eqn{p} and a discrete uniform variable \eqn{S} with probabilities \eqn{q_1,q_2,...,q_j} 
#' verifying \eqn{q_1+q_2+...+q_j=1-p}.
#'
#' @return Point and confidence estimates of the sensitive characteristics using the Eriksson model. The transformed variable is also reported, if required.
#'  
#' @references Eriksson, S.A. (1973). 
#' \emph{A new model for randomized response.}
#' International Statistical Review 41, 40-43.
#' 
#' @seealso \code{\link{ErikssonData}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Quantitative Eriksson Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=53376
#' data(ErikssonData)
#' dat=with(ErikssonData,data.frame(z,Pi))
#' p=0.5
#' mu=mean(c(0,1,3,5,8))
#' sigma=sqrt(4/5*var(c(0,1,3,5,8)))
#' cl=0.95
#' Eriksson(dat$z,p,mu,sigma,dat$Pi,"mean",cl,N)
#' 
#' @export
Eriksson=function(z,p,mu,sigma,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
 
  out1=Quantitative(z,p,0,1-p,c(0,0,mu),c(0,0,sigma))
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Eriksson",Model="Quantitative",Type=type,Param=list(c("p"=p,"mu"=mu,"sigma"=sigma)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}