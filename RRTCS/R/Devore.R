#' @name Devore
#' @aliases Devore
#' @title Devore model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Devore model. 
#' The function can also return the transformed variable. 
#' The Devore model was proposed by Devore in 1977.
#' 
#' @usage Devore(z,p,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p proportion of cards bearing the mark \eqn{A}
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level 
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' In the Devore model, the randomized response device presents to the sampled person labelled \eqn{i} a box containing a large number of identical cards with a 
#' proportion \eqn{p,(0<p<1)} bearing the mark \eqn{A} and the rest marked \eqn{B} (an innocuous attribute). The response solicited denoted by \eqn{z_i} takes the value 
#' \eqn{y_i} if \eqn{i} bears \eqn{A} and the card drawn is marked \eqn{A}. Otherwise \eqn{z_i} takes the value 1.
#' 
#' The transformed variable is \eqn{r_i=\frac{z_i-(1-p)}{p}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.
#'     
#' @return Point and confidence estimates of the sensitive characteristics using the Devore model. The transformed variable is also reported, if required.
#' 
#' @references Devore, J.L. (1977).
#' \emph{A note on the randomized response technique.}
#' Communications in Statistics Theory and Methods 6: 1525-1529.
#' 
#' @seealso \code{\link{DevoreData}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative Devore Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' data(DevoreData)
#' dat=with(DevoreData,data.frame(z,Pi))
#' p=0.7
#' cl=0.95
#' Devore(dat$z,p,dat$Pi,"total",cl)
#' 
#' @export
Devore=function(z,p,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,4,p,1)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Devore",Type=type,Model="Qualitative",Param=list(c("p"=p)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}