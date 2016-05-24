#' @name Kuk
#' @aliases Kuk
#' @title Kuk model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence through the Kuk model. 
#' The function can also return the transformed variable. 
#' The Kuk model was proposed by Kuk in 1990.
#' 
#' @usage Kuk(z,p1,p2,k,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p1 proportion of red cards in the first box 
#' @param p2 proportion of red cards in the second box  
#' @param k total number of cards drawn
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL
#' 
#' @details 
#' In the Kuk randomized response technique, the sampled person \eqn{i} is offered two boxes. Each box contains cards that are identical exception colour, either red 
#' or white, in sufficiently large numbers with proportions \eqn{p_1} and \eqn{1-p_1} in the first and \eqn{p_2} and \eqn{1-p_2}, in the second (\eqn{p_1\neq p_2}). 
#' The person sampled is requested to use the first box, if his/her trait is \eqn{A} and the second box if his/her trait is \eqn{A^c} and to make \eqn{k} independent 
#' draws of cards, with replacement each time. The person is asked to reports \eqn{z_i=f_i}, the number of times a red card is drawn.
#' 
#' The transformed variable is \eqn{r_i=\frac{f_i/k-p_2}{p_1-p_2}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=br_i+c}, 
#' where \eqn{b=\frac{1-p_1-p_2}{k(p_1-p_2)}} and \eqn{c=\frac{p_2(1-p_2)}{k(p_1-p_2)^2}}.
#'     
#' @return Point and confidence estimates of the sensitive characteristics using the Kuk model. The transformed variable is also reported, if required.
#'  
#' @references Kuk, A.Y.C. (1990). 
#' \emph{Asking sensitive questions indirectly.}
#'  Biometrika, 77, 436-438.
#'  
#' @seealso \code{\link{KukData}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative Kuk Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=802
#' data(KukData)
#' dat=with(KukData,data.frame(z,Pi))
#' p1=0.6
#' p2=0.2
#' k=25
#' cl=0.95
#' Kuk(dat$z,p1,p2,k,dat$Pi,"mean",cl,N)
#' 
#' @export
Kuk=function(z,p1,p2,k,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,7,p1,1,p2,k=k)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Kuk",Model="Qualitative",Type=type,Param=list(c("p1"=p1,"p2"=p2,"k"=k)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}