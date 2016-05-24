#' Compute the marginal MLE of phi
#'
#' Compute the marginal MLE of the ratio of two Poisson rates in a two-sample Poisson rate problem with misclassified data given fallible and infallible datasets.
#' 
#' @param data the vector of counts of the fallible data (z11, z12, z21, z22) followed by the infallible data (m011, m012, m021, m022, y01, y02)
#' @param N1 the opportunity size of group 1 for the fallible data
#' @param N2 the opportunity size of group 2 for the fallible data
#' @param N01 the opportunity size of group 1 for the infallible data
#' @param N02 the opportunity size of group 2 for the infallible data
#' @param l the lower end of the range of possible phi's (for optim)
#' @param u the upper end of the range of possible phi's (for optim)
#' @param out "par" or "all" (for the output of optim)
#' @return a named vector containing the marginal mle of phi
#' @export margMLE
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
#' margMLE(data, N1, N2, N01, N02)
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
#' margMLE(data, N1, N2, N01, N02)
#'
#'
#'
#'
#' }
#'
margMLE <- function(data, N1, N2, N01, N02, l = 1e-3, u = 1e3, out = c("par", "all")){
  
  out <- match.arg(out)
  
  phiHat <- fullMLE(un(data), N1, N2, N01, N02)[1]
  
  z11 <- data[1]; z12 <- data[2]; z21 <- data[3]; z22 <- data[4];
  m011 <- data[5]; m012 <- data[6]; m021 <- data[7]; m022 <- data[8]; 
  y01 <- data[9]; y02 <- data[10]
  
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
  
  if(o$convergence != 0) 
    warning('optim failed to converge; try out = "all" for details')
  
  if(out == "all") return(o)
  
  o$par
}
