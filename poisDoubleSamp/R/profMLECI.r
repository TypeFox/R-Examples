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
#' @export profMLECI
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
profMLECI <- function(data, N1, N2, N01, N02, conf.level = .95, l = 1e-3, u = 1e3, tol = 1e-10){
  
  # # define the function that, when 0, gives the boundaries of the CI  
  # ciFunction <- function(phi){
  #  mle <- un(fullMLE(data, N1, N2, N01, N02))
  #  
  #  unNormedLogLikelihood(
  #    data, N1, N2, N01, N02, mle[1], mle[2], mle[3], mle[4], mle[5], mle[6]
  #  ) - (
  #    unNormedLogProfileLikelihood(data, N1, N2, N01, N02, phi, tol) +
  #    qchisq(conf.level, df = 1)/2  
  #  )
  # }  
  
  ciFunction <- function(phi){
    mle <- un(fullMLE(data, N1, N2, N01, N02))
   
    logLikelihood(
      data, N1, N2, N01, N02, mle[1], mle[2], mle[3], mle[4], mle[5], mle[6]
    ) - (
      logProfileLikelihood(data, N1, N2, N01, N02, phi, tol) +
      qchisq(conf.level, df = 1)/2  
    )
  }  
  
  
  # compute limits and return
  phiHat <- un(fullMLE(data, N1, N2, N01, N02)["phi"])
  
  lower <- uniroot(ciFunction, lower = l, upper = phiHat, extendInt = "downX")$root
  
  upper <- uniroot(ciFunction, lower = phiHat, upper = u, extendInt = "upX")$root
  
  c(l = lower, u = upper)    
}


















eStep <- function(data, phiAndEta){
  
  data <- un(data)
   z11 <- data[1]; z12 <- data[2]; z21 <- data[3]; z22 <- data[4];
  m011 <- data[5]; m012 <- data[6]; m021 <- data[7]; m022 <- data[8]; 
   y01 <- data[9]; y02 <- data[10]
  
   phi <- un(phiAndEta["phi"])
  la12 <- un(phiAndEta["la12"]) 
  la21 <- un(phiAndEta["la21"])
  la22 <- un(phiAndEta["la22"])
   th1 <- un(phiAndEta["th1"])
   th2 <- un(phiAndEta["th2"])
  
  c(
    y1St = z12 * (phi*la21*th1) / (phi*la21*th1 + la12),
    y2St = z22 * (la21*th2) / (la21*th2 + la22)
  )    
}






mStep <- function(data, N1, N2, N01, N02, phi, completion){
  
  data <- un(data)
  
  z11 <- data[1]; z12 <- data[2]; z21 <- data[3]; z22 <- data[4];
  m011 <- data[5]; m012 <- data[6]; m021 <- data[7]; m022 <- data[8]; 
  y01 <- data[9]; y02 <- data[10]  
  
  y1St <- un(completion["y1St"]); y2St <- un(completion["y2St"])
  
  c(
    la12 = (z12 - y1St + m012) / (N1 + N01),
    la21 = (z11 + y1St + m011 + z21 + y2St + m021) / ((N1 + N01)*phi + (N2 + N02)),
    la22 = (z22 - y2St + m022) / (N2 + N02),
    th1 = (y1St + y01) / (z11 + y1St + m011),
    th2 = (y2St + y02) / (z21 + y2St + m021)
  )
}  








rmle <- function(data, N1, N2, N01, N02, phi, tol = 1e-10){        
  
  eta <- c(la12 = 1, la21 = 1, la22 = 1, th1 = .5, th2 = .5)
  
  # first step
  lastEta <- eta
  completion <- eStep(data, c(phi = phi, eta))
  eta <- mStep(data, N1, N2, N01, N02, phi, completion)
  
  # em algorihthm
  while(max(abs(eta - lastEta)) > tol){
    lastEta <- eta
    completion <- eStep(data, c(phi = phi, eta))
    eta <- mStep(data, N1, N2, N01, N02, phi, completion)
  }
  
  # return
  c(phi = phi, eta)
}

# rmle(data, N1, N2, N01, N02, 4)
# rmle(data, N1, N2, N01, N02, 5)
# fullMLE(data, N1, N2, N01, N02)
# rmle(data, N1, N2, N01, N02, 1.9129834)  # great!













unNormedLogLikelihood <- function(data, N1, N2, N01, N02, phi, la12, la21, la22, th1, th2){
  
  z11 <- data[1]; z12 <- data[2]; z21 <- data[3]; z22 <- data[4];
  m011 <- data[5]; m012 <- data[6]; m021 <- data[7]; m022 <- data[8]; 
  y01 <- data[9]; y02 <- data[10]  
  
  (m011 + m012)*log(N01) + (m021 + m022)*log(N02) +
    (z11 + z12)*log(N1) + (z21 + z22)*log(N2) +
    (m011 + z11)*log(phi) + (m011 + m021 + z11 + z21)*log(la21) +
    m012*log(la12) + m022*log(la22) + 
    y01*log(th1) + (m011 - y01 + z11)*log(1-th1) + 
    y02*log(th2) + (m021 - y02 + z21)*log(1-th2) + 
    z12*log(la12 + phi*la21*th1) + 
    z22*log(la22 + la21*th2) - 
    N01*phi*la21 - N01*la12 - N02*la21 - N02*la22 -
    N1*phi*la21*(1-th1) - N2*la21*(1-th2) -
    N1*(la12 + phi*la21*th1) - 
    N2*(la22 + la21*th2)
  
}







logLikelihood <- function(data, N1, N2, N01, N02, phi, la12, la21, la22, th1, th2){
  
  z11 <- data[1]; z12 <- data[2]; z21 <- data[3]; z22 <- data[4];
  m011 <- data[5]; m012 <- data[6]; m021 <- data[7]; m022 <- data[8]; 
  y01 <- data[9]; y02 <- data[10]  
  
  -(la12*(N01 + N1)) - (la21 + la22)*(N02 + N2) - la21*(N01 + N1)*phi + 
    m012*log(la12) + (m021 + z21)*log(la21) + m022*log(la22) + (m011 + m012)*log(N01) + 
    (m021 + m022)*log(N02) + (z11 + z12)*log(N1) + (z21 + z22)*log(N2) + 
    (m011 + z11)*(log(la21) + log(phi)) + (m011 - y01 + z11)*log(1 - th1) + 
    y01*log(th1) + z12*log(la12 + la21*phi*th1) + (m021 - y02 + z21)*log(1 - th2) + 
    y02*log(th2) + z22*log(la22 + la21*th2) + 
    lchoose(m011, y01) + lchoose(m021, y02) - 
    lfactorial(m011) - lfactorial(m012) - lfactorial(m021) - lfactorial(m022) - 
    lfactorial(z11) - lfactorial(z12) - lfactorial(z21) - lfactorial(z22)
}








unNormedLogProfileLikelihood <- function(data, N1, N2, N01, N02, phi, tol = 1e-10){

  phiAndEta <- rmle(data, N1, N2, N01, N02, phi, tol)
  
  la12 <- un(phiAndEta["la12"])
  la21 <- un(phiAndEta["la21"])
  la22 <- un(phiAndEta["la22"])
  th1 <- un(phiAndEta["th1"])
  th2 <- un(phiAndEta["th2"])
  
  unNormedLogLikelihood(data, N1, N2, N01, N02, phi, la12, la21, la22, th1, th2)
}
# fullMLE(data, N1, N2, N01, N02)
# unNormedLogLikelihood(data, N1, N2, N01, N02, 1.912983, 0.75, 2.784615, 2.523077, 0.422383, 0.226519)
# unNormedLogProfileLikelihood(data, N1, N2, N01, N02, 1.9129834)





logProfileLikelihood <- function(data, N1, N2, N01, N02, phi, tol = 1e-10){
  
  phiAndEta <- rmle(data, N1, N2, N01, N02, phi, tol)
  
  la12 <- un(phiAndEta["la12"])
  la21 <- un(phiAndEta["la21"])
  la22 <- un(phiAndEta["la22"])
  th1 <- un(phiAndEta["th1"])
  th2 <- un(phiAndEta["th2"])
  
  logLikelihood(data, N1, N2, N01, N02, phi, la12, la21, la22, th1, th2)
}
# fullMLE(data, N1, N2, N01, N02)
# logLikelihood(data, N1, N2, N01, N02, 1.912983, 0.75, 2.784615, 2.523077, 0.422383, 0.226519)
# logProfileLikelihood(data, N1, N2, N01, N02, 1.9129834)






