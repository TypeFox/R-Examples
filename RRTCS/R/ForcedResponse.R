#' @name ForcedResponse
#' @aliases ForcedResponse
#' @title Forced-Response model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Forced-Response model. 
#' The function can also return the transformed variable. 
#' The Forced-Response model was proposed by Boruch in 1972.
#' 
#' @usage ForcedResponse(z,p1,p2,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p1 proportion of cards marked "Yes"
#' @param p2 proportion of cards marked "No" 
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL
#' 
#' @details 
#' In the Forced-Response scheme, the sampled person \eqn{i} is offered a box with cards: some are marked "Yes" with a proportion \eqn{p_1}, some are marked "No" with a 
#' proportion \eqn{p_2} and the rest are marked "Genuine", in the remaining proportion \eqn{p_3=1-p_1-p_2}, where \eqn{0<p_1,p_2<1,p_1\neq p_2,p_1+p_2<1}. The person is 
#' requested to randomly draw one of them, to observe the mark on the card, and to respond
#' \deqn{z_i=\left \{\begin{array}{lccc}
#' 1 & \textrm{if the card is type "Yes"}\\
#' 0 & \textrm{if the card is type "No"}\\
#' y_i & \textrm{if the card is type "Genuine"}
#' \end{array}
#' \right .}
#' The transformed variable is \eqn{r_i=\frac{z_i-p_1}{1-p_1-p_2}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.  
#'
#' @return Point and confidence estimates of the sensitive characteristics using the Forced-Response model. The transformed variable is also reported, if required.
#' 
#' @references Boruch, R.F. (1972). 
#' \emph{Relations among statistical methods for assuring confidentiality of social research data.}
#'  Social Science Research, 1, 403-414.
#'   
#' @seealso \code{\link{ForcedResponseData}}
#' @seealso \code{\link{ForcedResponseDataSt}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative ForcedResponse Boruch Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' data(ForcedResponseData)
#' dat=with(ForcedResponseData,data.frame(z,Pi))
#' p1=0.2
#' p2=0.2
#' cl=0.95
#' ForcedResponse(dat$z,p1,p2,dat$Pi,"total",cl)
#' 
#' #Forced Response with strata
#' data(ForcedResponseDataSt)
#' dat=with(ForcedResponseDataSt,data.frame(ST,z,Pi))
#' p1=0.2
#' p2=0.2
#' cl=0.95
#' ForcedResponse(dat$z,p1,p2,dat$Pi,"total",cl)
#' 
#' @export 
ForcedResponse=function(z,p1,p2,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,6,p1,1,p2)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Forced Response",Model="Qualitative",Type=type,Param=list(c("p1"=p1,"p2"=p2)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}