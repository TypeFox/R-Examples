bayes.influence=function(theta,data)
{
y=data[,1]; n=data[,2]
N=length(y)
summary=quantile(theta[,2],c(.05,.5,.95))
summary.obs=array(0,c(N,3))
K=exp(theta[,2])
eta=exp(theta[,1])/(1+exp(theta[,1]))
m=length(K)

for (i in 1:N)
{
  weight=exp(lbeta(K*eta,K*(1-eta))-lbeta(K*eta+y[i],K*(1-eta)+n[i]-y[i]))
  probs=weight/sum(weight)
  indices=sample(1:m,size=m,prob=probs,replace=TRUE)
  theta.s=theta[indices,]
  summary.obs[i,]=quantile(theta.s[,2],c(.05,.5,.95))
}
return(list(summary=summary,summary.obs=summary.obs))
}
