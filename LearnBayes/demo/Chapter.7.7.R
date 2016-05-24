#########################################################
# Section 7.7 Simulating from the Posterior
#########################################################

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
 fitgibbs$accept

 windows()
 mycontour(poissgamexch, c(0, 8, -7.3, -6.6), datapar,
     xlab="log alpha",ylab="log mu")
 points(fitgibbs$par[, 1], fitgibbs$par[, 2])

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 plot(density(fitgibbs$par[, 1], bw = 0.2))

 alpha = exp(fitgibbs$par[, 1])
 mu = exp(fitgibbs$par[, 2])
 lam1 = rgamma(1000, y[1] + alpha, e[1] + alpha/mu)

 alpha = exp(fitgibbs$par[, 1])
 mu = exp(fitgibbs$par[, 2])

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 plot(log(e), y/e, pch = as.character(y))
 for (i in 1:94) {
     lami = rgamma(1000, y[i] + alpha, e[i] + alpha/mu)
     probint = quantile(lami, c(0.05, 0.95))
     lines(log(e[i]) * c(1, 1), probint)
 }


