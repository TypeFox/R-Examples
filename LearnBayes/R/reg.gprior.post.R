reg.gprior.post=function(theta,dataprior)
{
y=dataprior$data$y; X=dataprior$data$X
c0=dataprior$prior$c0; beta0=dataprior$prior$b0
beta=theta[-length(theta)]; sigma=exp(theta[length(theta)])

loglike=sum(dnorm(y,mean=X%*%as.vector(beta),sd=sigma,log=TRUE))
logprior=dmnorm(beta,mean=beta0,varcov=c0*sigma^2*solve(t(X)%*%X),log=TRUE)
return(loglike+logprior)
}