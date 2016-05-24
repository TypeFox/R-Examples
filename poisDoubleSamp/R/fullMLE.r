#' Compute the full MLEs
#'
#' Compute the MLEs of a two-sample Poisson rate problem with misclassified data given fallible and infallible datasets.
#' 
#' These are the closed-form expressions for the MLEs.
#' 
#' @param data the vector of counts of the fallible data (z11, z12, z21, z22) followed by the infallible data (m011, m012, m021, m022, y01, y02)
#' @param N1 the opportunity size of group 1 for the fallible data
#' @param N2 the opportunity size of group 2 for the fallible data
#' @param N01 the opportunity size of group 1 for the infallible data
#' @param N02 the opportunity size of group 2 for the infallible data
#' @return a named vector containing the mles of each of the parameters (phi, la12, la21, la22, th1, and th2)
#' @export fullMLE
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
#' fullMLE(data, N1, N2, N01, N02)
#' 
#' 
#' 
#' # big example :
#' z11 <- 477; z12 <- 1025; N1 <- 16186;
#' z21 <- 255; z22 <- 1450; N2 <- 18811;
#' m011 <- 38;  m012 <- 90; y01 <- 15; N01 <- 1500; 
#' m021 <- 41; m022 <- 200; y02 <-  9; N02 <- 2500;
#' data <- c(z11, z12, z21, z22, m011, m012, m021, m022, y01, y02)
#' 
#' fullMLE(data, N1, N2, N01, N02)
#'
#'
#' }
#'
fullMLE <- function(data, N1, N2, N01, N02){
  
  data <- un(data)
  z11 <- data[1]; z12 <- data[2]; z21 <- data[3]; z22 <- data[4];
  m011 <- data[5]; m012 <- data[6]; m021 <- data[7]; m022 <- data[8]; 
  y01 <- data[9]; y02 <- data[10]
  
  phiHat <- exp(
    log(N2 + N02) + log(m022 + y02) + 
      log(m011*m012 + m011*y01 + m012*z11 + y01*z11 + y01*z12) -
      log(N1 + N02) - log(m012 + y01) - 
      log(m021*m022 + m021*y02 + m022*z21 + y02*z21 + y02*z22)
  )
  
  la21Hat <- exp(
    log(m021*m022 + m021*y02 + m022*z21 + y02*z21 + y02*z22) -
      log(N2 + N02) - log(m022 + y02)
  )
  
  la12Hat <- (m012^2 + m012*y01 + m012*z12) / (N1 + N01) / (m012 + y01)
  
  la22Hat <- (m022^2 + m022*y02 + m022*z22) / (N2 + N02) / (m022 + y02)
  
  th1Hat <- 1 - (m012 + y01) * (m011 - y01 + z11) / 
    (m011*m012 + m011*y01 + m012*z11 + y01*z11 + y01*z12)
  
  th2Hat <- 1 - (m022 + y02) * (m021 - y02 + z21) / 
    (m021*m022 + m021*y02 + m022*z21 + y02*z21 + y02*z22)
  
  c(
    phi = phiHat, la12 = la12Hat, la21 = la21Hat, la22 = la22Hat,
    th1 = th1Hat, th2 = th2Hat
  )   
}
