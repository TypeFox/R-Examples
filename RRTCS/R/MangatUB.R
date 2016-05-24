#' @name MangatUB
#' @aliases MangatUB
#' @title Mangat-UB model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Mangat model (Mangat, 1992)
#' when the proportion of people bearing the innocuous attribute is unknown. 
#' The function can also return the transformed variable. 
#' The Mangat-UB model can be seen in Chaudhuri (2011, page 53).
#' 
#' @usage MangatUB(I,J,p1,p2,t,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param I first vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param J second vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p1 proportion of marked cards with the sensitive attribute in the second box
#' @param p2 proportion of marked cards with the sensitive attribute in the third box
#' @param t probability of response to the sensitive questions without using random response in the first box
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' In Mangat's extended scheme, three boxes containing cards are presented to the sampled person, labelled \eqn{i}. The first box contains cards marked "True" and "RR" 
#' in proportions \eqn{t} and \eqn{1-t}, the second one contains \eqn{A} and \eqn{B}-marked cards in proportions \eqn{p_1} and \eqn{(1-p_1),(0<p_1<1)} and the third box 
#' contains \eqn{A} and \eqn{B}-marked cards in proportions \eqn{p_2} and \eqn{1-p_2,(0<p_2<1),p_1\neq p_2}. The subject is requested to draw a card from the first box. 
#' The sample respondent \eqn{i} is then instructed to tell the truth, using "the first box and if necessary also the second box" and next, independently, to give a second 
#' truthful response also using "the first box and if necessary, the third box." Let \eqn{I_i} represent the first response and \eqn{J_i} the second response for 
#' respondent \eqn{i}.
#' 
#' The transformed variable is \eqn{r_i=\frac{(1-p_2)I_i-(1-p_1)J_i}{p_1-p_2}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.
#' 
#' @return Point and confidence estimates of the sensitive characteristics using the Mangat-UB model. The transformed variable is also reported, if required.
#' 
#' @references Chaudhuri, A. (2011). 
#' \emph{Randomized response and indirect questioning techniques in surveys.}
#' Boca Raton: Chapman and Hall, CRC Press.  
#' 
#' @references Mangat, N.S. (1992). 
#' \emph{Two stage randomized response sampling procedure using unrelated question.}
#' Journal of the Indian Society of Agricultural Statistics, 44, 82-87.
#'  
#' @seealso \code{\link{Mangat}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative Mangat Estimation Variance Transformed_variable Confidence_interval
#' 
#' @export
MangatUB=function(I,J,p1,p2,t,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(I,13,p1,2,p2,t=t,z2=J)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Mangat unknown B",Model="Qualitative",Type=type,Param=list(c("p1"=p1,"p2"=p2,"t"=t)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}