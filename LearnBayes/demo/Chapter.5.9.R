#############################################
# Section 5.9 Importance Sampling
#############################################

library(LearnBayes)

data(cancermortality)

fit=laplace(betabinexch,c(-7,6),cancermortality)

betabinexch.cond=function (log.K, data)
{
eta = exp(-6.818793)/(1 + exp(-6.818793))
K = exp(log.K)
y = data[, 1]; n = data[, 2]; N = length(y)
logf=0*log.K
for (j in 1:length(y))
logf = logf + lbeta(K * eta + y[j], K * (1 -
eta) + n[j] - y[j]) - lbeta(K * eta, K * (1 - eta))
val = logf + log.K - 2 * log(1 + K)
return(exp(val-max(val)))
}


I=integrate(betabinexch.cond,2,16,cancermortality)
par(mfrow=c(2,2))
curve(betabinexch.cond(x,cancermortality)/I$value,from=3,to=16,
ylab="Density", xlab="log K",lwd=3, main="Densities")
curve(dnorm(x,8,2),add=TRUE)
legend("topright",legend=c("Exact","Normal"),lwd=c(3,1))
curve(betabinexch.cond(x,cancermortality)/I$value/
  dnorm(x,8,2),from=3,to=16, ylab="Weight",xlab="log K",
  main="Weight = g/p")

curve(betabinexch.cond(x,cancermortality)/I$value,from=3,to=16,
   ylab="Density", xlab="log K",lwd=3, main="Densities")
curve(1/2*dt(x-8,df=2),add=TRUE)
legend("topright",legend=c("Exact","T(2)"),lwd=c(3,1))
curve(betabinexch.cond(x,cancermortality)/I$value/
        (1/2*dt(x-8,df=2)),from=3,to=16, ylab="Weight",xlab="log K",
        main="Weight = g/p")

tpar=list(m=fit$mode,var=2*fit$var,df=4)
 myfunc=function(theta)
   return(theta[2])
 s=impsampling(betabinexch,tpar,myfunc,10000,cancermortality)
 cbind(s$est,s$se)

