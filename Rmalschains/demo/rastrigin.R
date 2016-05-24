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

res.rastrigin1 <- malschains(rastrigin, lower=rep(-1.0, 30), upper=rep(1.0, 30), maxEvals=50000, 
    control=malschains.control(effort=0.8, alpha=0.3, popsize=20, istep=100, ls="simplex"))

# We supply an initialpop of one individual. The other ones will be initialized randomly
res.rastrigin2 <- malschains(rastrigin, lower=rep(-1.0, 30), upper=rep(1.0, 30), maxEvals=50000, 
    initialpop = rep(0.1, 30), control=malschains.control(popsize=50, istep=300, ls="cmaes"))

res.rastrigin3 <- malschains(rastrigin, lower=rep(-1.0, 30), upper=rep(1.0, 30), maxEvals=50000, 
    control=malschains.control(effort=0.8, alpha=0.3, popsize=50, istep=100, ls="sw"))

res.rastrigin1
res.rastrigin2
res.rastrigin3

