#' Compute the Wald confidence interval
#'
#' Compute the Wald confidence interval of a two-sample Poisson rate with misclassified data given fallible and infallible datasets.
#' 
#' @param data the vector of counts of the fallible data (z11, z12, z21, z22) followed by the infallible data (m011, m012, m021, m022, y01, y02)
#' @param N1 the opportunity size of group 1 for the fallible data
#' @param N2 the opportunity size of group 2 for the fallible data
#' @param N01 the opportunity size of group 1 for the infallible data
#' @param N02 the opportunity size of group 2 for the infallible data
#' @param conf.level confidence level of the interval
#' @return a named vector containing the lower and upper bounds of the confidence interval
#' @export waldCI
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
#' }
#'
waldCI <- function(data, N1, N2, N01, N02, conf.level = .95){
  
  mle <- unname(fullMLE(data, N1, N2, N01, N02))
  se <- sePhi(mle, data, N1, N2, N01, N02)
  alpha <- 1-conf.level
  
  c(
    l = mle[1] - qnorm(1-alpha/2) * se,
    u = mle[1] + qnorm(1-alpha/2) * se
  )
}
























sePhi <- function(params, data, N1, N2, N01, N02){
  
  i <- infoMatrix(params, data, N1, N2, N01, N02)
  
  1 / sqrt(as.numeric(i[1,1] - i[2:6,1] %*% solve(i[2:6, 2:6]) %*% i[1,2:6]))
}












infoMatrix <- function(params, data, N1, N2, N01, N02){

  params <- un(params)
  data <- un(data)
  
  phi <- params[1];
  la12 <- params[2]; la21 <- params[3]; la22 <- params[4];
  th1 <- params[5]; th2 <- params[6]
  
  z11 <- data[1]; z12 <- data[2]; z21 <- data[3]; z22 <- data[4];
  m011 <- data[5]; m012 <- data[6]; m021 <- data[7]; m022 <- data[8]; 
  y01 <- data[9]; y02 <- data[10]
  
  info <- matrix(0, nrow = 6, ncol = 6)
  
  info[1,1] <- (m011 + z11)/phi^2 + (la21^2*th1^2*z12)/(la12 + la21*phi*th1)^2
  info[1,2] <- (la21*th1*z12)/(la12 + la21*phi*th1)^2
  info[1,3] <- N01 + N1 - (la12*th1*z12)/(la12 + la21*phi*th1)^2
  info[1,4] <- 0
  info[1,5] <- -((la12*la21*z12)/(la12 + la21*phi*th1)^2)
  info[1,6] <- 0
  info[2,1] <- (la21*th1*z12)/(la12 + la21*phi*th1)^2
  info[2,2] <- m012/la12^2 + z12/(la12 + la21*phi*th1)^2
  info[2,3] <- (phi*th1*z12)/(la12 + la21*phi*th1)^2
  info[2,4] <- 0
  info[2,5] <- (la21*phi*z12)/(la12 + la21*phi*th1)^2
  info[2,6] <- 0
  info[3,1] <- N01 + N1 - (la12*th1*z12)/(la12 + la21*phi*th1)^2
  info[3,2] <- (phi*th1*z12)/(la12 + la21*phi*th1)^2
  info[3,3] <- (phi^2*th1^2*z12)/(la12 + la21*phi*th1)^2 + (m011 + m021 + z11 + z21)/la21^2 + (th2^2*z22)/(la22 + la21*th2)^2
  info[3,4] <- (th2*z22)/(la22 + la21*th2)^2
  info[3,5] <- -((la12*phi*z12)/(la12 + la21*phi*th1)^2)
  info[3,6] <- -((la22*z22)/(la22 + la21*th2)^2)
  info[4,1] <- 0
  info[4,2] <- 0
  info[4,3] <- (th2*z22)/(la22 + la21*th2)^2
  info[4,4] <- m022/la22^2 + z22/(la22 + la21*th2)^2
  info[4,5] <- 0
  info[4,6] <- (la21*z22)/(la22 + la21*th2)^2
  info[5,1] <- -((la12*la21*z12)/(la12 + la21*phi*th1)^2)
  info[5,2] <- (la21*phi*z12)/(la12 + la21*phi*th1)^2
  info[5,3] <- -((la12*phi*z12)/(la12 + la21*phi*th1)^2)
  info[5,4] <- 0
  info[5,5] <- y01/th1^2 + (m011 - y01 + z11)/(-1 + th1)^2 + (la21^2*phi^2*z12)/(la12 + la21*phi*th1)^2
  info[5,6] <- 0
  info[6,1] <- 0
  info[6,2] <- 0
  info[6,3] <- -((la22*z22)/(la22 + la21*th2)^2)
  info[6,4] <- (la21*z22)/(la22 + la21*th2)^2
  info[6,5] <- 0
  info[6,6] <- y02/th2^2 + (m021 - y02 + z21)/(-1 + th2)^2 + (la21^2*z22)/(la22 + la21*th2)^2    
  
  info
}


