#' @name SinghJoarder
#' @aliases SinghJoarder
#' @title Singh-Joarder model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Singh-Joarder model. 
#' The function can also return the transformed variable. 
#' The Singh-Joarder model was proposed by Singh and Joarder in 1997.
#' 
#' @usage SinghJoarder(z,p,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p proportion of marked cards with the sensitive question 
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level 
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' The basics of the Singh-Joarder scheme are similar to Warner's randomized response device, with the following difference. If a person labelled \eqn{i} bears 
#' \eqn{A^c} he/she is told to say so if so guided by a card drawn from a box of \eqn{A} and \eqn{A^c} marked cards in proportions \eqn{p} and \eqn{(1-p),(0<p<1)}.
#' However, if he/she bears \eqn{A} and is directed by the card to admit it, he/she is told to postpone the reporting based on the first draw of the card from 
#' the box but to report on the basis of a second draw. Therefore,
#' \deqn{z_i=\left \{\begin{array}{lcc}
#' 1 & \textrm{if person } i \textrm{ responds "Yes"}\\
#' 0 & \textrm{if person } i \textrm{ responds "No"}
#' \end{array}
#' \right .}   
#' The transformed variable is \eqn{r_i=\frac{z_i-(1-p)}{(2p-1)+p(1-p)}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.
#'
#' @return Point and confidence estimates of the sensitive characteristics using the Singh-Joarder model. The transformed variable is also reported, if required.
#'  
#' @references Singh, S., Joarder, A.H. (1997).
#' \emph{Unknown repeated trials in randomized response sampling.}
#' Journal of the Indian Statistical Association, 30, 109-122.
#' 
#' @seealso \code{\link{SinghJoarderData}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative SinghJoarder Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=802
#' data(SinghJoarderData)
#' dat=with(SinghJoarderData,data.frame(z,Pi))
#' p=0.6
#' cl=0.95
#' SinghJoarder(dat$z,p,dat$Pi,"mean",cl,N)
#' 
#' @export
SinghJoarder=function(z,p,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,9,p,1)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Singh and Joarder",Model="Qualitative",Type=type,Param=list(c("p"=p)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out) 
}