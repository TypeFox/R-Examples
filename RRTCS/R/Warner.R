#' @name Warner
#' @aliases Warner
#' @title Warner model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Warner model. 
#' The function can also return the transformed variable. 
#' The Warner model was proposed by Warner in 1965.
#' 
#' @usage Warner(z,p,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p proportion of marked cards with the sensitive attribute  
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level 
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' Warner's randomized response device works as follows. A sampled person labelled \eqn{i} is offered a box of a considerable number of identical cards with a proportion 
#' \eqn{p,(0<p<1,p\neq 0.5)} of them marked \eqn{A} and the rest marked \eqn{A^c}. The person is requested, randomly, to draw one of them, to observe the mark on the card, 
#' and to give the response
#' \deqn{z_i=\left\{\begin{array}{lcc}
#'  1 & \textrm{if card type "matches" the trait } A \textrm{ or } A^c \\
#'  0 & \textrm{if a "no match" results }
#'  \end{array}
#'  \right.}
#'  The randomized response is given by \eqn{r_i=\frac{z_i-(1-p)}{2p-1}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.
#'     
#' @return Point and confidence estimates of the sensitive characteristics using the Warner model. The transformed variable is also reported, if required.
#' 
#' @references Warner, S.L. (1965). 
#' \emph{Randomized Response: a survey technique for eliminating evasive answer bias.}
#'  Journal of the American Statistical Association 60, 63-69.
#'  
#' @seealso \code{\link{WarnerData}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative Warner Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=802
#' data(WarnerData)
#' dat=with(WarnerData,data.frame(z,Pi))
#' p=0.7
#' cl=0.95
#' Warner(dat$z,p,dat$Pi,"total",cl)
#' 
#' @export
Warner=function(z,p,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,1,p,1)
  out2=Estimator(out1,pi,type,cl,N,pij)

  out=c(Call=call,out1,out2,Name="Warner",Model="Qualitative",Type=type,Param=list(c("p"=p)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}