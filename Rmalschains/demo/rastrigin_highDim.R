library(Rmalschains)

rastrigin <- function(x) {
  
  dimension <- length(x)
  
  res <- 0.0
  for (i in 1:dimension) {
    res <- res + (x[i]*x[i] - 10.0*cos(2.0*pi*x[i]) + 10.0)
  }

  #print(paste("fitness:", res, sep=""))
  #res <- res - 330
  res 
}

res.rastrigin.highDim <- malschains(rastrigin, lower=rep(-1.0, 1000), upper=rep(1.0, 1000), maxEvals=50000, 
    control=malschains.control(effort=0.8, alpha=0.3, popsize=20, istep=100, ls="ssw"))

res.rastrigin.highDim
