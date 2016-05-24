#' @name ChaudhuriChristofides
#' @aliases ChaudhuriChristofides
#' @title Chaudhuri-Christofides model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Chaudhuri-Christofides model. 
#' The function can also return the transformed variable. 
#' The Chaudhuri-Christofides model can be seen in Chaudhuri and Christofides (2013, page 97).
#' 
#' @usage ChaudhuriChristofides(z,mu,sigma,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param mu vector with the means of the scramble variables
#' @param sigma vector with the standard deviations of the scramble variables
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL 
#' 
#' @details 
#' The randomized response given by the person \eqn{i} is \eqn{z_i=y_iS_1+S_2} where \eqn{S_1,S_2} are scramble variables, whose mean \eqn{\mu} and 
#' standard deviation \eqn{\sigma} are known. 
#'
#' @return Point and confidence estimates of the sensitive characteristics using the Chaudhuri-Christofides model. The transformed variable is also reported, if required.
#' 
#' @references Chaudhuri, A., and Christofides, T.C. (2013)
#' \emph{Indirect Questioning in Sample Surveys.}
#' Springer-Verlag Berlin Heidelberg.
#' 
#' @seealso \code{\link{ChaudhuriChristofidesData}}
#' @seealso \code{\link{ChaudhuriChristofidesDatapij}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Quantitative ChaudhuriChristofides Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=417
#' data(ChaudhuriChristofidesData)
#' dat=with(ChaudhuriChristofidesData,data.frame(z,Pi))
#' mu=c(6,6) 
#' sigma=sqrt(c(10,10))
#' cl=0.95
#' data(ChaudhuriChristofidesDatapij)
#' ChaudhuriChristofides(dat$z,mu,sigma,dat$Pi,"mean",cl,pij=ChaudhuriChristofidesDatapij)
#' 
#' @export
ChaudhuriChristofides=function(z,mu,sigma,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
 
  out1=Quantitative(z,0,1,0,c(mu,0),c(sigma,0))
  out2=Estimator(out1,pi,type,cl,N,pij)
  
  out=c(Call=call,out1,out2,Name="Chaudhuri and Christofides",Model="Quantitative",Type=type,Param=list(c("mu"=mu,"sigma"=sigma)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}