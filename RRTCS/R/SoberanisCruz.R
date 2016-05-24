#' @name SoberanisCruz
#' @aliases SoberanisCruz
#' @title SoberanisCruz model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the SoberanisCruz model. 
#' The function can also return the transformed variable. 
#' The SoberanisCruz model was proposed by Soberanis Cruz et al. in 2008.
#' 
#' @usage SoberanisCruz(z,p,alpha,pi,type=c("total","mean"),cl,N=NULL,pij=NULL)
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param p proportion of marked cards with the sensitive question  
#' @param alpha proportion of people with the innocuous attribute
#' @param pi vector of the first-order inclusion probabilites
#' @param type the estimator type: total or mean 
#' @param cl confidence leve
#' @param N size of the population. By default it is NULL
#' @param pij matrix of the second-order inclusion probabilities. By default it is NULL
#' 
#' @details 
#' The SoberanisCruz model considers the introduction of an innocuous variable correlated with the sensitive variable. 
#' This variable does not affect individual sensitivity, and maintains reliability.
#' The sampling procedure is the same as in the Horvitz model.
#' 
#' @return Point and confidence estimates of the sensitive characteristics using the SoberanisCruz model. The transformed variable is also reported, if required.
#' 
#' @references Soberanis Cruz, V., Ramírez Valverde, G., Pérez Elizalde, S., González Cossio, F. (2008).
#' \emph{Muestreo de respuestas aleatorizadas en poblaciones finitas: Un enfoque unificador.} 
#'  Agrociencia Vol. 42 Núm. 5 537-549.
#'  
#' @seealso \code{\link{SoberanisCruzData}}
#' @seealso \code{\link{Horvitz}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Qualitative SoberanisCruz Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' data(SoberanisCruzData)
#' dat=with(SoberanisCruzData,data.frame(z,Pi))
#' p=0.7
#' alpha=0.5 
#' cl=0.90
#' SoberanisCruz(dat$z,p,alpha,dat$Pi,"total",cl)
#' 
#' @export
SoberanisCruz=function(z,p,alpha,pi,type=c("total","mean"),cl,N=NULL,pij=NULL){
  call=match.call()
  
  out1=Qualitative(z,3,p,1,alpha=alpha)
  out2=Estimator(out1,pi,type,cl,N,pij)
 
  out=c(Call=call,out1,out2,Name="Soberanis Cruz",Model="Qualitative",Type=type,Param=list(c("p"=p,"alpha"=alpha)),ConfidenceLevel=cl)
  class(out)="RRT"
  return(out)
}