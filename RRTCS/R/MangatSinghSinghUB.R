#' @name MangatSinghSinghUB
#' @aliases MangatSinghSinghUB
#' @title Mangat-Singh-Singh-UB model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Mangat-Singh-Singh model (Mangat el al., 1992)
#' when the proportion of people bearing the innocuous attribute is unknown. 
#' The function can also return the transformed variable. 
#' The Mangat-Singh-Singh-UB model can be seen in Chauduri (2011, page 54).
#' 
#' @usage MangatSinghSinghUB(I,J,p1,p2,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param I first vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param J second vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p1 proportion of marked cards with the sensitive attribute in the first box
#' @param p2 proportion of marked cards with the sensitive attribute in the second box
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' A person labelled \eqn{i} who is chosen, is instructed to say "yes" if he/she bears \eqn{A}, and if not, to randomly take a card from a box containing
#' cards marked \eqn{A,B} in proportions \eqn{p_1} and \eqn{(1-p_1),(0<p_1<1)}; they are then told to report the value \eqn{x_i} if a \eqn{B}-type card is chosen and he/she bears \eqn{B}; 
#' otherwise he/she is told to report "No". This entire exercise is to be repeated independently with the second box with \eqn{A} and \eqn{B}-marked cards in proportions \eqn{p_2} 
#' and \eqn{(1-p_2),(0<p_2<1,p_2\neq p_1)}. Let \eqn{I_i} the first response and \eqn{J_i} the second response for the respondent \eqn{i}.
#' 
#' The transformed variable is \eqn{r_i=\frac{(1-p_2)I_i-(1-p_1)J_i}{p_1-p_2}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}. 
#'
#' @return Point and confidence estimates of the sensitive characteristics using the Mangat-Singh-Singh-UB model. The transformed variable is also reported, if required.
#' 
#' @references Chaudhuri, A. (2011). 
#' \emph{Randomized response and indirect questioning techniques in surveys.}
#' Boca Raton: Chapman and Hall, CRC Press.  
#'  
#' @references Mangat, N.S., Singh, R., Singh, S. (1992). 
#' \emph{An improved unrelated question randomized response strategy.}
#' Calcutta Statistical Association Bulletin, 42, 277-281.
#'  
#' @seealso \code{\link{MangatSinghSinghUBData}}
#' @seealso \code{\link{MangatSinghSingh}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative MangatSinghSingh Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=802
#' data(MangatSinghSinghUBData)
#' dat=with(MangatSinghSinghUBData,data.frame(I,J,Pi))
#' p1=0.6
#' p2=0.8
#' cl=0.95
#' MangatSinghSinghUB(dat$I,dat$J,p1,p2,dat$Pi,"mean",cl,N)
#' 
#' @export
MangatSinghSinghUB=function(I,J,p1,p2,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(I,14,p1,2,p2,z2=J)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Mangat, Singh and Singh unknown B",Model="Qualitative",Type=type,Param=list(c("p1"=p1,"p2"=p2)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}