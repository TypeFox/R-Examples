#' Compute the profile MLE CI of phi
#'
#' Compute the profile MLE confidence interval of the ratio of two Poisson rates in a two-sample Poisson rate problem with misclassified data given fallible and infallible datasets.  This uses a C++ implemention of the EM algorithm.
#' 
#' @param data the vector of counts of the fallible data (z11, z12, z21, z22) followed by the infallible data (m011, m012, m021, m022, y01, y02)
#' @param N1 the opportunity size of group 1 for the fallible data
#' @param N2 the opportunity size of group 2 for the fallible data
#' @param N01 the opportunity size of group 1 for the infallible data
#' @param N02 the opportunity size of group 2 for the infallible data
#' @param conf.level confidence level of the interval
#' @param l the lower end of the range of possible phi's (for optim)
#' @param u the upper end of the range of possible phi's (for optim)
#' @param tol tolerance used in the EM algorithm to declare convergence
#' @return a named vector containing the marginal mle of phi
#' @export approxMargMLECI
#' @examples
#' \dontrun{
#'
#' # small example
#' z11 <- 34; z12 <- 35; N1 <- 10; 
#' z21 <- 22; z22 <- 31; N2 <- 10;
#' m011 <- 9; m012 <- 1; y01 <- 3; N01 <- 3;
#' m021 <- 8; m022 <- 8; y02 <- 2; N02 <- 3;
#' data <- c(z11, z12, z21, z22, m011, m012, m021, m022, y01, y02)
#' 
#' waldCI(data, N1, N2, N01, N02) 
#' margMLECI(data, N1, N2, N01, N02)
#' profMLECI(data, N1, N2, N01, N02)
#' approxMargMLECI(data, N1, N2, N01, N02)
#' 
#' 
#' # big example :
#' z11 <- 477; z12 <- 1025; N1 <- 16186;
#' z21 <- 255; z22 <- 1450; N2 <- 18811;
#' m011 <- 38;  m012 <- 90; y01 <- 15; N01 <- 1500; 
#' m021 <- 41; m022 <- 200; y02 <-  9; N02 <- 2500;
#' data <- c(z11, z12, z21, z22, m011, m012, m021, m022, y01, y02)
#' 
#' waldCI(data, N1, N2, N01, N02) 
#' margMLECI(data, N1, N2, N01, N02)
#' profMLECI(data, N1, N2, N01, N02)
#' approxMargMLECI(data, N1, N2, N01, N02)
#'
#'
#'
#' }
#'
approxMargMLECI <- function(data, N1, N2, N01, N02, conf.level = .95, l = 1e-3, u = 1e3, tol = 1e-10){

  approxMargPhiHat <- approxMargMLE(data, N1, N2, N01, N02, l, u, tol = tol)  
  
  # define the function that, when 0, gives the boundaries of the CI  
  #zeroFunc <- function(phi){
  #  unNormedApproxMargLogLikelihood(data, N1, N2, N01, N02, phi, tol) - 
  #    unNormedApproxMargLogLikelihood(data, N1, N2, N01, N02, approxMargPhiHat, tol) +
  #    .5 * qchisq(conf.level, df = 1)/2  
  #}
  
  ciFunction <- function(phi){
    approxMargLogLikelihood(data, N1, N2, N01, N02, phi, tol) - 
      approxMargLogLikelihood(data, N1, N2, N01, N02, approxMargPhiHat, tol) +
      .5 * qchisq(conf.level, df = 1)/2  
  }
  
  # compute limits and return
  lower <- uniroot(ciFunction, lower = l, upper = approxMargPhiHat, 
    extendInt = "upX")$root
  
  upper <- uniroot(ciFunction, lower = approxMargPhiHat, upper = u, 
    extendInt = "downX")$root
  
  c(l = lower, u = upper)    
}
