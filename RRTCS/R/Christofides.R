#' @name Christofides
#' @aliases Christofides
#' @title Christofides model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Christofides model. 
#' The function can also return the transformed variable. 
#' The Christofides model was proposed by Christofides in 2003.
#' 
#' @usage Christofides(z,mm,pm,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param mm vector with the marks of the cards
#' @param pm vector with the probabilities of previous marks 
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level 
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' In the Christofides randomized response technique, a sampled person \eqn{i} is given a box with identical cards, each bearing a separate mark as 
#' \eqn{1,\dots,k,\dots m} with \eqn{m\geq 2} but in known proportions \eqn{p_1,\dots,p_k,\dots p_m} with \eqn{0<p_k< 1} for \eqn{k=1,\dots,m} and 
#' \eqn{\sum_{k=1}^{m}p_k=1}. The person sampled is requested to draw one of the cards and respond
#' \deqn{z_i=\left \{\begin{array}{lcc}
#' k & \textrm{if a card marked } k \textrm{ is drawn and the person bears } A^c\\
#' m-k+1 & \textrm{if a card marked } k \textrm{ is drawn but the person bears } A
#' \end{array}
#' \right .}   
#' The transformed variable is \eqn{r_i=\frac{z_i-\mu}{m+1-2\mu}} where \eqn{\mu=\sum_{k=1}^{m}kp_k} and the estimated variance is
#' \eqn{\widehat{V}_R(r_i)=\frac{V_R(k)}{(m+1-2\mu)^2}}, where \eqn{V_R(k)=\sum_{k=1}^{m}k^2p_k-\mu^2}.
#'
#' @return Point and confidence estimates of the sensitive characteristics using the Christofides model. The transformed variable is also reported, if required.
#'  
#' @references Christofides, T.C. (2003). 
#' \emph{A generalized randomized response technique.}
#'  Metrika, 57, 195-200.
#'  
#' @seealso \code{\link{ChristofidesData}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative Christofides Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=802
#' data(ChristofidesData)
#' dat=with(ChristofidesData,data.frame(z,Pi))
#' mm=c(1,2,3,4,5)
#' pm=c(0.1,0.2,0.3,0.2,0.2)
#' cl=0.95
#' Christofides(dat$z,mm,pm,dat$Pi,"mean",cl,N)
#' 
#' @export
Christofides=function(z,mm,pm,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,8,nrr=1,mm=mm,pm=pm)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Christofides",Model="Qualitative",Type=type,Param=list(c("mm"=mm,"pm"=pm)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}