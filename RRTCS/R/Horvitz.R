#' @name Horvitz
#' @aliases Horvitz
#' @title Horvitz model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Horvitz model. 
#' The function can also return the transformed variable. 
#' The Horvitz model was proposed by Horvitz et al. (1967) and by Greenberg et al. (1969).
#' 
#' @usage Horvitz(z,p,alpha,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p proportion of marked cards with the sensitive question 
#' @param alpha proportion of people with the innocuous attribute
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level 
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL
#' 
#' @details 
#' In the Horvitz model, the randomized response device presents to the sampled person labelled \eqn{i} a box containing a large number of identical cards, with a 
#' proportion \eqn{p,(0 <p < 1)} bearing the mark \eqn{A} and the rest marked \eqn{B} (an innocuous attribute whose population proportion \eqn{\alpha} is known). 
#' The response solicited denoted by \eqn{z_i} takes the value \eqn{y_i} if \eqn{i} bears \eqn{A} and the card drawn is marked \eqn{A} or if \eqn{i} bears 
#' \eqn{B} and the card drawn is marked \eqn{B}. Otherwise \eqn{z_i} takes the value 0.
#' 
#' The transformed variable is \eqn{r_i=\frac{z_i-(1-p)\alpha}{p}} and the estimated variance is \eqn{\widehat{V}_R(r_i)=r_i(r_i-1)}.
#' 
#' @return Point and confidence estimates of the sensitive characteristics using the Horvitz model. The transformed variable is also reported, if required.
#'  
#' @references Greenberg, B.G., Abul-Ela, A.L., Simmons, W.R., Horvitz, D.G. (1969).
#' \emph{The unrelated question RR model: Theoretical framework.}
#' Journal of the American Statistical Association, 64, 520-539.
#' 
#' @references Horvitz, D.G., Shah, B.V., Simmons, W.R. (1967).
#' \emph{The unrelated question RR model.}
#'  Proceedings of the Social Statistics Section of the American Statistical Association. 65-72. Alexandria, VA: ASA.
#' 
#' @seealso \code{\link{HorvitzData}}
#' @seealso \code{\link{HorvitzDataStCl}}
#' @seealso \code{\link{HorvitzDataRealSurvey}}
#' @seealso \code{\link{HorvitzUB}}
#' @seealso \code{\link{SoberanisCruz}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative Horvitz Greenberg Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=10777
#' data(HorvitzData)
#' dat=with(HorvitzData,data.frame(z,Pi))
#' p=0.5
#' alpha=0.6666667
#' cl=0.95
#' Horvitz(dat$z,p,alpha,dat$Pi,"mean",cl,N) 
#' 
#' #Horvitz real survey
#' N=10777
#' n=710
#' data(HorvitzDataRealSurvey)
#' p=0.5
#' alpha=1/12
#' pi=rep(n/N,n)
#' cl=0.95
#' Horvitz(HorvitzDataRealSurvey$sex,p,alpha,pi,"mean",cl,N)
#' 
#' @export
Horvitz=function(z,p,alpha,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,2,p,1,alpha=alpha)
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Horvitz",Model="Qualitative",Type=type,Param=list(c("p"=p,"alpha"=alpha)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}