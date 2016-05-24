StephensTestRad <- function(x, mu=0, alpha) {
  n <- length(x)
  Muobs <- MeanCircularRad(x)
  Robs <- RhoCircularRad(x)*n
  Cobs <- Robs*cos(Muobs-mu) 
  Rcrit <- Ralpha(x=Cobs, n=n, alpha=alpha)
  test <- TRUE
  if (Robs > Rcrit) test <- FALSE
  result <- c(test, Robs, Cobs, Rcrit)
  return(result)
}
