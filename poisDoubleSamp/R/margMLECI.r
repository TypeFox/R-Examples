#' Compute the marginal MLE confidence interval for the phi
#'
#' Compute the marginal MLE confidence interval of the ratio of two Poisson rates in a two-sample Poisson rate problem with misclassified data given fallible and infallible datasets.
#' 
#' @param data the vector of counts of the fallible data (z11, z12, z21, z22) followed by the infallible data (m011, m012, m021, m022, y01, y02)
#' @param N1 the opportunity size of group 1 for the fallible data
#' @param N2 the opportunity size of group 2 for the fallible data
#' @param N01 the opportunity size of group 1 for the infallible data
#' @param N02 the opportunity size of group 2 for the infallible data
#' @param conf.level confidence level of the interval
#' @param l the lower end of the range of possible phi's (for optim)
#' @param u the upper end of the range of possible phi's (for optim)
#' @return a named vector containing the lower and upper bounds of the confidence interval
#' @export margMLECI
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
margMLECI <- function(data, N1, N2, N01, N02, conf.level = .95, l = 1e-10, u = 1e10){
  
  z11 <- data[1]; z12 <- data[2]; z21 <- data[3]; z22 <- data[4];
  m011 <- data[5]; m012 <- data[6]; m021 <- data[7]; m022 <- data[8]; 
  y01 <- data[9]; y02 <- data[10];
  
  phiHat <- fullMLE(data, N1, N2, N01, N02)[1]
  
  maxLogTerm <- unNormedMargLogLikeMaxLogTermCpp(
    phiHat, z11, z12, z21, z22, m011, m012, m021, m022, 
    y01, y02, N1, N2, N01, N02
  )
  
  Lmarg <- function(x) unNormedMargLogLikeCpp(x, 
    maxLogTerm, z11, z12, z21, z22, m011, m012, m021, m022, 
    y01, y02, N1, N2, N01, N02
  )  
  
  #o <- optim(phiHat, fn = Lmarg, control = list(fnscale = -1),
  #      method = "Brent", lower = l, upper = u)
  # found to be sensitive for some problems
  
  o <- optim(phiHat, fn = Lmarg, control = list(fnscale = -1), method = "BFGS")
  
  margPhiHat <- o$par
  lAtMargPhiHat <- o$value
  
  unNormLikeMinusVal <- function(phi)
    unNormedMargLogLikeCpp(phi, 
      maxLogTerm, z11, z12, z21, z22, m011, m012, m021, m022, 
      y01, y02, N1, N2, N01, N02
    ) - (lAtMargPhiHat - qchisq(conf.level, df = 1)/2)
  
  
  
  lower <- suppressWarnings(uniroot(unNormLikeMinusVal, 
    lower = l, upper = margPhiHat, extendInt = "upX")$root)
  
  upper <- suppressWarnings(uniroot(unNormLikeMinusVal, 
    lower = margPhiHat, upper = u, extendInt = "downX")$root)
  
  c(l = lower, u = upper)
}
