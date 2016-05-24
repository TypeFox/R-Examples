#' @name Mangat
#' @aliases Mangat
#' @title Mangat model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Mangat model. 
#' The function can also return the transformed variable. 
#' The Mangat model was proposed by Mangat in 1992.
#' 
#' @usage Mangat(z,p,alpha,t,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p proportion of marked cards with the sensitive attribute in the second box
#' @param alpha proportion of people with the innocuous attribute
#' @param t proportion of marked cards with "True" in the first box
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' In Mangat's method, there are two boxes, the first containing cards marked "True" and "RR" in proportions \eqn{t} and \eqn{(1-t),(0<t<1)}. A person drawing a 
#' "True" marked card is asked to tell the truth about bearing \eqn{A} or \eqn{A^c}. A person drawing and “RR” marked card is then asked to apply Horvitz’s device 
#' by drawing a card from a second box with cards marked \eqn{A} and \eqn{B} in proportions \eqn{p} and \eqn{(1-p)}. If an \eqn{A} marked card is now drawn the 
#' truthful response will be about bearing the sensitive attribute \eqn{A} and otherwise about \eqn{B}. The true proportion of people bearing \eqn{A} is to be 
#' estimated but \eqn{\alpha}, the proportion of people bearing the innocuous trait \eqn{B} unrelated to \eqn{A}, is assumed to be known. The observed variable is
#' \deqn{z_i=\left \{\begin{array}{lcc}
#' y_i & \textrm{if a card marked "True" is drawn from the first box}\\
#' I_i & \textrm{if a card marked "RR" is drawn}
#' \end{array}
#' \right .}
#' where
#' \deqn{I_i=\left \{\begin{array}{lcc}
#' 1 & \textrm{if the type of card drawn from the second box matches trait } A \textrm{ or } B\\
#' 0 & \textrm{if the type of card drawn from the second box does not match trait } A \textrm{ or } B.
#' \end{array}
#' \right .}
#' 
#' The transformed variable is \eqn{r_i=\frac{z_i-(1-t)(1-p)\alpha}{t+(1-t)p}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.
#'     
#' @return Point and confidence estimates of the sensitive characteristics using the Mangat model. The transformed variable is also reported, if required.
#' 
#' @references Mangat, N.S. (1992). 
#' \emph{Two stage randomized response sampling procedure using unrelated question.}
#' Journal of the Indian Society of Agricultural Statistics, 44, 82-87.
#'  
#' @seealso \code{\link{MangatUB}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative Mangat Estimation Variance Transformed_variable Confidence_interval
#' 
#' @export
Mangat=function(z,p,alpha,t,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,10,p,1,alpha=alpha,t=t)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Mangat",Model="Qualitative",Type=type,Param=list(c("p"=p,"alpha"=alpha,"t"=t)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}