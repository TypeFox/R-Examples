#'  Hellinger distance between two MCMC chains using a specified grid size.
#'  
#'  This function determines the Hellinger distance between two MCMC chains via kernel density estimates.
#'  
#'  @param b1 vector of first MCMC chain.
#'  @param b2 vector of second MCMC chain.
#'  @param n2 is the number of divisions to run for the kernel density estimator
#'  @return res1 The Hellinger distance between the kernel density estimates for b1 and b2.
#'  @note The chains need to be the same length.
HDistSize <- function(b1,b2,n2){
  #n1 <- nrow(as.matrix(b1))
  b2c <- c(b1,b2)
  b2min <- min(b2c)
  b2max <- max(b2c)
  P1 <- density(b1,from=b2min,to=b2max,n=n2)
  Q1 <- density(b2,from=b2min,to=b2max,n=n2)
  Pdiff1 <- P1$y
  Qdiff1 <- Q1$y
  step1 <- P1$x[2]-P1$x[1]
  diver1 <- (sqrt(Pdiff1)-sqrt(Qdiff1))^2*step1
  res1 <- sqrt(sum(diver1)/2)
  return(res1)
}

#' MCMC.one is an mcmc object resulting from the following code:
#'
#' @name MCMC.one
#' @docType data
#' @author Edward L. Boone \email{elboone@@vcu.edu}
#' @keywords data
#' @examples
#' \dontrun{
#'  library(dismo); library(MCMCpack)
#'  data(Anguilla_train)
#'  b0mean0 <- 0
#'  b0precision <- (1/5)^2
#'  mcmclen = 50000
#'  burn=200000
#'  MCMC.one <- MCMClogit(Angaus ~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
#'                   data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-1, 
#'                   b0=b0mean0, B0=b0precision)
#'  }
NULL

#' MCMC.two is an mcmc object resulting from the following code:
#'
#' @name MCMC.two
#' @docType data
#' @author Edward L. Boone \email{elboone@@vcu.edu}
#' @keywords data
#' @examples
#' \dontrun{
#'  library(dismo); library(MCMCpack)
#'  data(Anguilla_train)
#'  b0mean0 <- 0
#'  b0precision <- (1/5)^2
#'  mcmclen = 50000
#'  burn=200000
#'  MCMC.two <- MCMClogit(Angaus ~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
#'                   data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-1, 
#'                   b0=b0mean0, B0=b0precision)
#'  }
NULL

#' MCMC.three is an mcmc object resulting from the following code:
#'
#' @name MCMC.three
#' @docType data
#' @author Edward L. Boone \email{elboone@@vcu.edu}
#' @keywords data
#' @examples
#' \dontrun{
#'  library(dismo); library(MCMCpack)
#'  data(Anguilla_train)
#'  b0mean0 <- 0
#'  b0precision <- (1/5)^2
#'  mcmclen = 50000
#'  burn=200000
#'  MCMC.three <- MCMClogit(Angaus ~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
#'                   data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-1, 
#'                   b0=b0mean0, B0=b0precision)
#'  }
NULL

#' MCMC.one.mean0 is an mcmc object resulting from the following code:
#'
#' @name MCMC.one.mean0
#' @docType data
#' @author Edward L. Boone \email{elboone@@vcu.edu}
#' @keywords data
#' @examples
#' \dontrun{
#'  library(dismo); library(MCMCpack)
#'  data(Anguilla_train)
#'  b0mean0 <- 0
#'  b0precision <- (1/5)^2
#'  mcmclen = 50000
#'  burn=200000
#'  MCMC.one.mean0 <- MCMClogit(Angaus ~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
#'                   data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-1, 
#'                   b0=b0mean0, B0=b0precision)
#'  }
NULL

#' MCMC.one.mean1 is an mcmc object resulting from the following code:
#'
#' @name MCMC.one.mean1
#' @docType data
#' @author Edward L. Boone \email{elboone@@vcu.edu}
#' @keywords data
#' @examples
#' \dontrun{
#'  library(dismo); library(MCMCpack)
#'  data(Anguilla_train)
#'  b0mean1 <- 1
#'  b0precision <- (1/5)^2
#'  mcmclen = 50000
#'  burn=200000
#'  MCMC.one.mean1 <- MCMClogit(Angaus ~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
#'                   data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-1, 
#'                   b0=b0mean1, B0=b0precision)
#'  }
NULL