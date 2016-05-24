robustt=function(y,v,m)
{
rigamma=function(n,a,b)
{
# simulates n values from a Inverse Gamma
# distribution with shape a and rate b
# density x^(-a-1) exp(b/x)

return(1/rgamma(n,shape=a,rate=b))
}
n=length(y)
mu=mean(y); sig2=sd(y)^2; lam=array(1,c(n,1))
M=array(0,c(m,1)); S2=M; LAM=array(0,c(m,n))

for (i in 1:m)
{
  lam=rgamma(n,shape=(v+1)/2,rate=v/2+(y-mu)^2/2/sig2)
  mu=rnorm(1,mean=sum(y*lam)/sum(lam),sd=sqrt(sig2/sum(lam)))
  sig2=rigamma(1,n/2,sum(lam*(y-mu)^2)/2)
  M[i]=mu; S2[i]=sig2; LAM[i,]=lam
}
par=list(mu=M,s2=S2,lam=LAM)
return(par)	  			 
}
