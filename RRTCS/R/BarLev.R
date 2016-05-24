#' @name BarLev
#' @aliases BarLev
#' @title BarLev model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the BarLev model.
#' The function can also return the transformed variable. 
#' The BarLev model was proposed by Bar-Lev et al. in 2004.
#' 
#' @usage BarLev(z,p,mu,sigma,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p probability of direct response 
#' @param mu mean of the scramble variable \eqn{S}
#' @param sigma standard deviation of the scramble variable \eqn{S}
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' The randomized  response given by the person \eqn{i} is 
#' \deqn{z_i=\left\{\begin{array}{lcc}
#' y_i & \textrm{with probability } p\\
#' y_iS & \textrm{with probability } 1-p\\
#' \end{array}
#' \right.}
#' where \eqn{S} is a scramble variable, whose mean \eqn{\mu} and standard deviation \eqn{\sigma} are known. 
#'  
#' @return Point and confidence estimates of the sensitive characteristics using the BarLev model. The transformed variable is also reported, if required.
#'
#' @references Bar-Lev S.K., Bobovitch, E., Boukai, B. (2004). 
#' \emph{A note on randomized response models for quantitative data.}
#' Metrika, 60, 255-260.
#' 
#' @seealso \code{\link{BarLevData}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Quantitative BarLev Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' data(BarLevData)
#' dat=with(BarLevData,data.frame(z,Pi))
#' p=0.6
#' mu=1
#' sigma=1
#' cl=0.95
#' BarLev(dat$z,p,mu,sigma,dat$Pi,"total",cl)
#' 
#' @export
BarLev=function(z,p,mu,sigma,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Quantitative(z,p,1-p,0,c(mu,0,0),c(sigma,0,0))
  out2=Estimator(out1,pi,type,cl,N,pij)
 
  out=c(Call=call,out1,out2,Name="Bar Lev",Model="Quantitative",Type=type,Param=list(c("p"=p,"mu"=mu,"sigma"=sigma)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}