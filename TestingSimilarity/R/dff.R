################################################################################
#' Implementation of absolute difference function
#'
#' Function calculating the absolute difference of two dose response models: \deqn{dff(d,\theta_1,\theta_2)=|m_1(d,\theta_1)-m_2(d,\theta_2)|}
#' 
#' @export
#' @param d real-valued argument to the function (dose variable)
#' @param alpha,beta model parameters (real vectors)
#' @param m1,m2 model types. Built-in models are "linlog",  "linear",  "quadratic",  "emax",  "exponential",  "sigEmax",  "betaMod" and "logistic" 
#' @return Response value for the absolute difference of two models.
################################################################################  
dff <- function(d,alpha,beta,m1,m2){
  abs(m1(d,alpha)-m2(d,beta))
}