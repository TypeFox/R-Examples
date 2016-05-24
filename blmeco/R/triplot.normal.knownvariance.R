triplot.normal.knownvariance <-
function(theta.data, variance.known, n, prior.theta, prior.variance, 
         legend=TRUE, ylim=c(0, max(yposterior)), legend.bty="n"){
  
x <- rnorm(n, mean=theta.data, sd=sqrt(variance.known))
xlim <- theta.data+c(-1, 1)*3*sqrt(variance.known)
steplength <- (xlim[2]-xlim[1])/500
xx <- seq(xlim[1], xlim[2], by=steplength)
yprior <- dnorm(xx, mean=prior.theta, sd=sqrt(prior.variance))

likelihood.function <- function(x, theta, variance.known) prod(1/(sqrt(2*pi*sqrt(variance.known)))*exp(-1/(2*variance.known)*(x-theta)^2))
likelihood <- numeric(length(xx))
for(i in 1:length(xx)) likelihood[i] <- likelihood.function(x, xx[i], variance.known=variance.known)

posterior.mean <- (prior.theta/prior.variance + n*mean(x)/variance.known)/(1/prior.variance + n/variance.known)
posterior.variance <- 1/(1/prior.variance + n/variance.known)
yposterior  <- dnorm(xx, posterior.mean, sqrt(posterior.variance))

linewidth <- 2
plot(xx, yprior, type="l", lty=3, las=1, ylab="density",
  xlab=expression(theta), cex.lab=1.4, cex.axis=1.1, lwd=linewidth, ylim=ylim)
lines(xx, yposterior, lty=1, lwd=linewidth)
lines(xx, likelihood/sum(likelihood)*1/steplength, lty=2, lwd=linewidth)
if(legend) legend(xlim[1], ylim[2], lwd=linewidth, lty=c(3,2,1), 
                  legend=c("prior", "likelihood", "posterior"), bty=legend.bty)
return(list(posterior.mean=posterior.mean, posterior.variance=posterior.variance, x=x))
}
