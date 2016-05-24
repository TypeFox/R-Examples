#' @name MangatSingh
#' @aliases MangatSingh
#' @title Mangat-Singh model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Mangat-Singh model. 
#' The function can also return the transformed variable. 
#' The Mangat-Singh model was proposed by Mangat and Singh in 1990.
#' 
#' @usage MangatSingh(z,p,t,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p proportion of marked cards with the sensitive attribute in the second box
#' @param t proportion of marked cards with "True" in the first box
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' In the Mangat-Singh model, the sampled person is offered two boxes of cards. In the first box a known proportion \eqn{t,(0<t<1)} of cards is marked "True" and the 
#' remaining ones are marked "RR". One card is to be drawn, observed and returned to the box. If the card drawn is marked "True", then the respondent should respond "Yes" 
#' if he/she belongs to the sensitive category, otherwise "No". If the card drawn is marked "RR", then the respondent must use the second box and draw a card from it. 
#' This second box contains a proportion \eqn{p,(0<p<1,p\neq 0.5)} of cards marked \eqn{A} and the remaining ones are marked \eqn{A^c}. If the card drawn from the second 
#' box matches his/her status as related to the stigmatizing characteristic, he/she must respond "Yes", otherwise "No".
#' The randomized response from a person labelled \eqn{i} is assumed to be:
#' \deqn{z_i=\left \{\begin{array}{lcc}
#' y_i & \textrm{if a card marked "True" is drawn from the first box}\\
#' I_i & \textrm{if a card marked "RR" is drawn}
#' \end{array}
#' \right .}  
#' \deqn{I_i=\left \{\begin{array}{lcc}
#' 1 & \textrm{if the "card type" } A \textrm{ or } A^c  \textrm{ "matches" the genuine trait } A \textrm{ or } A^c\\
#' 0 & \textrm{if a "mismatch" is observed}
#' \end{array}
#' \right .}
#' The transformed variable is \eqn{r_i=\frac{z_i-(1-t)(1-p)}{t+(1-t)(2p-1)}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.
#'  
#' @return Point and confidence estimates of the sensitive characteristics using the Mangat-Singh model. The transformed variable is also reported, if required.
#' 
#' @references Mangat, N.S., Singh, R. (1990). 
#' \emph{An alternative randomized response procedure.}
#'  Biometrika, 77, 439-442.
#' 
#' @seealso \code{\link{MangatSinghData}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative MangatSingh Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=802
#' data(MangatSinghData)
#' dat=with(MangatSinghData,data.frame(z,Pi))
#' p=0.7
#' t=0.55
#' cl=0.95
#' MangatSingh(dat$z,p,t,dat$Pi,"mean",cl,N)
#' 
#' @export
MangatSingh=function(z,p,t,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,5,p,1,t=t)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Mangat and Singh",Model="Qualitative",Type=type,Param=list(c("p"=p,"t"=t)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}