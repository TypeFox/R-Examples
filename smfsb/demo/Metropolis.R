# Metropolis.R

# Simple Metropolis sampler with different parameter settings

require(smfsb)

n=10000
alphas=c(0.1,1,100)

p=length(alphas)
op=par(mfrow=c(p,3))
for (alpha in alphas) {
    vec=metrop(n,alpha)
    plot(ts(vec),xlab="Iteration",ylab="Value",main=paste("Alpha =",alpha))
    acf(vec,lag.max=100,main="ACF")
    hist(vec,30,xlim=c(-4,4),xlab="Value",main="Empirical density",freq=FALSE)
}
par(op)


# eof


