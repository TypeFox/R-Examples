poisson.gamma.mix=function(probs,gammapar,data)
{
N=length(probs)
y=data$y; t=data$t; n=length(y)
post.gammapar=gammapar+outer(rep(1,N),c(sum(y),sum(t)))
L=post.gammapar[,1]/post.gammapar[,2]

loglike=0
for (j in 1:n)
  loglike=loglike+dpois(y[j],L*t[j],log=TRUE)

m.prob=exp(loglike+
  dgamma(L,shape=gammapar[,1],rate=gammapar[,2],log=TRUE) -
  dgamma(L,shape=post.gammapar[,1],rate=post.gammapar[,2],log=TRUE))

post.probs=probs*m.prob/sum(probs*m.prob)
return(list(probs=post.probs,gammapar=post.gammapar))
}