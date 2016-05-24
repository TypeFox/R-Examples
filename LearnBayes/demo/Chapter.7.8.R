##########################################################
# Section 7.8  Posterior Inferences
##########################################################

 library(LearnBayes)
 data(hearttransplants)
 attach(hearttransplants)

 datapar = list(data = hearttransplants, z0 = 0.53)
 start=c(2, -7)
 fit = laplace(poissgamexch, start, datapar)
 fit

 par(mfrow = c(1, 1))
 mycontour(poissgamexch, c(0, 8, -7.3, -6.6), datapar,
   xlab="log alpha",ylab="log mu")

S=readline(prompt="Type  <Return>   to continue : ")

 start = c(4, -7)
 fitgibbs = gibbs(poissgamexch, start, 1000, c(1,.15), datapar)

alpha = exp(fitgibbs$par[, 1])
 mu = exp(fitgibbs$par[, 2])

shrink=function(i) mean(alpha/(alpha + e[i] * mu))
 shrinkage=sapply(1:94, shrink)

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 plot(log(e), shrinkage)

 mrate=function(i) mean(rgamma(1000, y[i] + alpha, e[i] + alpha/mu))
 hospital=1:94
 meanrate=sapply(hospital,mrate)
 hospital[meanrate==min(meanrate)]

###########################################################

sim.lambda=function(i) rgamma(1000,y[i]+alpha,e[i]+alpha/mu)
LAM=sapply(1:94,sim.lambda)

compare.rates <- function(x) {
  nc <- NCOL(x)
  ij <- as.matrix(expand.grid(1:nc, 1:nc))
  m <- as.matrix(x[,ij[,1]] > x[,ij[,2]]) 
  matrix(colMeans(m), nc, nc, byrow = TRUE)
}

better=compare.rates(LAM)

better[1:24,85]
