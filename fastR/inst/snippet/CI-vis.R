# simulate 100 intervals and plot them.
results <- CIsim(n=20, samples=100, estimand=500, 
	rdist=rnorm, args=list(mean=500,sd=100),
	method=ci, method.args=list(sd=100))
coverPlot <- xYplot(Cbind(estimate,lower,upper) ~ sample, results,
	groups=cover, col=c('black','gray40'),cap=0,lwd=2,pch=16)
