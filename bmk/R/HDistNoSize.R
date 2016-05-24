#'  Hellinger distance between two MCMC chains using default grid in kernel density estimator.
#'  
#'  This function determines the Hellinger distance between two MCMC chains via kernel density estimates.
#'  
#'  @param b1 vector of first MCMC chain.
#'  @param b2 vector of second MCMC chain.
#'  @return The Hellinger distance between the kernel density estimates for b1 and b2.
#'  @note The chains need to be the same length.
HDistNoSize <- function(b1,b2){
  n1 <- nrow(as.matrix(b1))
  b2c <- c(b1,b2)
  b2min <- min(b2c)
  b2max <- max(b2c)
  P1 <- density(b1,from=b2min,to=b2max,n=n1)
  Q1 <- density(b2,from=b2min,to=b2max,n=n1)
  Pdiff1 <- P1$y
  Qdiff1 <- Q1$y
  step1 <- P1$x[2]-P1$x[1]
  diver1 <- (sqrt(Pdiff1)-sqrt(Qdiff1))^2*step1
  res1 <- sqrt(sum(diver1)/2)
  return(res1)
}
