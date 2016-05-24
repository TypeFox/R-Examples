plot.ssym <-
function(x, ...){
par(mfrow=c(2,2))

if(x$censored==FALSE) plot(x$z_es, x$weights, main = "Standardized individual-specific weights", cex=0.3, lwd=3, xlab="Ordinary residuals", ylab="Weight")
else plot(exp(x$z_es), x$cdfz, main = "Survival function of the error", cex=0.3, lwd=3, xlab="Multiplicative Error", ylab="Survival Function")

if(x$censored==FALSE) qqnorm(qnorm(x$cdfz), main="Overall goodness-of-fit statistic", cex=0.3, lwd=3, xlab="Quantiles of N(0,1)", ylab="Overall residuals")
else{surv0 <- survfit(Surv(x$z_es,1-x$event)~1)
ids <- ifelse(surv0$surv>0 & surv0$surv<1 & surv0$n.event>0,TRUE,FALSE)
plot(qnorm(1-surv0$surv[ids]),qnorm(x$cdf(surv0$time[ids])), main="Overall goodness-of-fit statistic", cex=0.3, lwd=3, xlab="Empirical Distribution", ylab="Expected Distribution")}
abline(0,1)

res.dev.mu <- sqrt(x$deviance.mu)*ifelse(x$z_es>=0,1,-1)   
ry <- c(min(res.dev.mu,-3.5),max(res.dev.mu,3.5))               
if(x$censored==FALSE) plot(x$mu.fitted,res.dev.mu, ylim=ry, cex=0.3, lwd=3, main="Median submodel", ylab="Deviance-type residual", xlab="Fitted values")                                 
else{plot(x$mu.fitted[x$event==0],res.dev.mu[x$event==0], xlim=range(x$mu.fitted), ylim=ry, cex=0.3, lwd=3, main="Median submodel", ylab="Deviance-type residual", xlab="Fitted values")
par(new=TRUE)
plot(x$mu.fitted[x$event==1],res.dev.mu[x$event==1], xlim=range(x$mu.fitted), ylim=ry, pch="+", main="Median submodel", ylab="Deviance-type residual", xlab="Fitted values")} 
abline(h=-3,lty=3)                                              
abline(h=+3,lty=3)                                              

res.dev.phi <- sqrt(x$deviance.phi)*ifelse(x$z_es>=0,1,-1) 
ry <- c(min(res.dev.phi,-3.5),max(res.dev.phi,3.5))             
if(x$censored==FALSE) plot(x$phi.fitted,res.dev.phi, ylim=ry, cex=0.3, lwd=3, main="Skewness/Dispersion submodel", ylab="Deviance-type residual", xlab="Fitted values")                                 
else{plot(x$phi.fitted[x$event==0],res.dev.phi[x$event==0], xlim=range(x$phi.fitted), ylim=ry, cex=0.3, lwd=3, main="Skewness/Dispersion submodel", ylab="Deviance-type residual", xlab="Fitted values")
par(new=TRUE)
plot(x$phi.fitted[x$event==1],res.dev.phi[x$event==1], xlim=range(x$phi.fitted), ylim=ry, pch="+", main="Skewness/Dispersion submodel", ylab="Deviance-type residual", xlab="Fitted values")} 
abline(h=-3,lty=3)                                              
abline(h=+3,lty=3)
}
