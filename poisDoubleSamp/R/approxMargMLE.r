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
#' @param tol tolerance parameter for the rmle EM algorithm
#' @return a named vector containing the marginal mle of phi
#' @export approxMargMLE
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
#' approxMargMLE(data, N1, N2, N01, N02)
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
#' # margMLE(data, N1, N2, N01, N02) # ~1 min
#' approxMargMLE(data, N1, N2, N01, N02)
#'
#'
#'
#' }
#'
approxMargMLE <- function(data, N1, N2, N01, N02, l = 0, u = 1e3, out = c("par", "all"), tol = 1e-10){
  
  out <- match.arg(out)
  
  phiHat <- fullMLE(data, N1, N2, N01, N02)[1]
  
  funcToBeMaximized <- function(phi) 
    unNormedApproxMargLogLikelihood(data, N1, N2, N01, N02, phi, tol)
  
  
  o <- optim(phiHat, fn = funcToBeMaximized, 
    control = list(fnscale = -1),
    method = "Brent", lower = l, upper = u)
  
  if(o$convergence != 0) 
    warning('optim failed to converge; try out = "all" for details')
  
  if(out == "all") return(o)
  
  o$par
}















logDetJeta <- function(data, N1, N2, N01, N02, phi, la12, la21, la22, th1, th2){

  ieta <- infoMatrix(
    c(phi, la12, la21, la22, th1, th2), 
    data, 
    N1, N2, N01, N02
  )[2:6,2:6]
  
  d <- determinant(ieta, logarithm = TRUE)
  
  as.numeric(d$sign * d$modulus)
}
#logDetJeta(data, N1, N2, N01, N02, phi, la12, la21, la22, th1, th2)









unNormedApproxMargLogLikelihood <- function(data, N1, N2, N01, N02, phi, tol = 1e-10){
  
  phiAndEta <- rmle(data, N1, N2, N01, N02, phi, tol)
  
  la12Hat <- un(phiAndEta["la12"])
  la21Hat <- un(phiAndEta["la21"])
  la22Hat <- un(phiAndEta["la22"])
  th1Hat <- un(phiAndEta["th1"])
  th2Hat <- un(phiAndEta["th2"])

  
  # we could have used the unNormedLogProfileLikelihood function, but then
  # we'd have to compute the rmle twice, so we do it this way
  unNormedLogLikelihood(data, N1, N2, N01, N02, 
                        phi, la12Hat, la21Hat, la22Hat, th1Hat, th2Hat) - 
    .5 * logDetJeta(data, N1, N2, N01, N02,
                    phi, la12Hat, la21Hat, la22Hat, th1Hat, th2Hat)
}







approxMargLogLikelihood <- function(data, N1, N2, N01, N02, phi, tol = 1e-10){
  
  phiAndEta <- rmle(data, N1, N2, N01, N02, phi, tol)
  
  la12Hat <- un(phiAndEta["la12"])
  la21Hat <- un(phiAndEta["la21"])
  la22Hat <- un(phiAndEta["la22"])
  th1Hat <- un(phiAndEta["th1"])
  th2Hat <- un(phiAndEta["th2"])
  
  
  # we could have used the logProfileLikelihood function, but then
  # we'd have to compute the rmle twice, so we do it this way
  5/2*log(2*pi) +
    logLikelihood(data, N1, N2, N01, N02, phi, la12Hat, la21Hat, la22Hat, th1Hat, th2Hat) - 
    .5 * logDetJeta(data, N1, N2, N01, N02, phi, la12Hat, la21Hat, la22Hat, th1Hat, th2Hat)
}
