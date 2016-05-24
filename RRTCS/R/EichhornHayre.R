#' @name EichhornHayre
#' @aliases EichhornHayre
#' @title Eichhorn-Hayre model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Eichhorn-Hayre model. 
#' The function can also return the transformed variable. 
#' The Eichhorn-Hayre model was proposed by Eichhorn and Hayre in 1983.
#' 
#' @usage EichhornHayre(z,mu,sigma,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param mu mean of the scramble variable \eqn{S}
#' @param sigma standard deviation of the scramble variable \eqn{S}
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level 
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' The randomized response given by the person labelled \eqn{i} is \eqn{z_i=y_iS} where \eqn{S} is a scramble variable whose distribution is assumed to be known.
#'
#' @return Point and confidence estimates of the sensitive characteristics using the Eichhorn-Hayre model. The transformed variable is also reported, if required.
#'   
#' @references  Eichhorn, B.H., Hayre, L.S. (1983). 
#' \emph{Scrambled randomized response methods for obtaining sensitive quantitative data.}
#' Journal of Statistical Planning and Inference, 7, 306-316.
#'  
#' @seealso \code{\link{EichhornHayreData}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Quantitative EichhornHayre Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' data(EichhornHayreData)
#' dat=with(EichhornHayreData,data.frame(z,Pi))
#' mu=1.111111
#' sigma=0.5414886
#' cl=0.95
#' #This line returns a warning showing why the variance estimation is not possible. 
#' #See ResamplingVariance for several alternatives. 
#' EichhornHayre(dat$z,mu,sigma,dat$Pi,"mean",cl)
#'  
#' @export
EichhornHayre=function(z,mu,sigma,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Quantitative(z,0,1,0,c(mu,0,0),c(sigma,0,0))
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Eichhorn and Hayre",Model="Quantitative",Type=type,Param=list(c("mu"=mu,"sigma"=sigma)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}