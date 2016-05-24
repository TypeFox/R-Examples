library(Rmalschains)

sphere <- function(x)
{
  
  sum <- 0.0
  for (i in 1:length(x)) {
    sum <- sum + abs(x[i]+i)  
  }
  
  #print(paste("x: ", x, sep=""))  
  #print(paste("sum: ", sum, sep=""))
  
  sum
}

res.sphere <- malschains(sphere, lower=rep(-1.0, 30), upper=rep(1.0, 30), maxEvals=50000, 
    control=malschains.control(effort=0.8, alpha=0.3, popsize=20, istep=100, ls="simplex"))

res.sphere

