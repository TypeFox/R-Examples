#' @name MangatSinghSingh
#' @aliases MangatSinghSingh
#' @title Mangat-Singh-Singh model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Mangat-Singh-Singh model. 
#' The function can also return the transformed variable. 
#' The Mangat-Singh-Singh model was proposed by Mangat, Singh and Singh in 1992.
#' 
#' @usage MangatSinghSingh(z,p,alpha,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p proportion of marked cards with the sensitive attribute in the box
#' @param alpha proportion of people with the innocuous attribute
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' In the Mangat-Singh-Singh scheme, a person labelled \eqn{i}, if sampled, is offered a box and told to answer "yes" if the person bears \eqn{A}. But if the person bears 
#' \eqn{A^c} then the person is to draw a card from the box with a proportion \eqn{p(0<p< 1)} of cards marked \eqn{A} and the rest marked \eqn{B}; if the person draws 
#' a card marked \eqn{B} he/she is told to say "yes" again if he/she actually bears \eqn{B}; in any other case, "no" is to be answered.
#' 
#' The transformed variable is \eqn{r_i=\frac{z_i-(1-p)\alpha}{1-(1-p)\alpha}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.
#'     
#' @return Point and confidence estimates of the sensitive characteristics using the Mangat-Singh-Singh model. The transformed variable is also reported, if required.
#' 
#' @references Mangat, N.S., Singh, R., Singh, S. (1992). 
#' \emph{An improved unrelated question randomized response strategy.}
#' Calcutta Statistical Association Bulletin, 42, 277-281.
#' 
#' @seealso \code{\link{MangatSinghSinghData}}
#' @seealso \code{\link{MangatSinghSinghUB}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative MangatSinghSingh Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' data(MangatSinghSinghData)
#' dat=with(MangatSinghSinghData,data.frame(z,Pi))
#' p=0.6
#' alpha=0.5
#' cl=0.95 
#' MangatSinghSingh(dat$z,p,alpha,dat$Pi,"total",cl)
#' 
#' @export
MangatSinghSingh=function(z,p,alpha,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,11,p,1,alpha=alpha)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Mangat,Singh and Singh",Model="Qualitative",Type=type,Param=list(c("p"=p,"alpha"=alpha)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}