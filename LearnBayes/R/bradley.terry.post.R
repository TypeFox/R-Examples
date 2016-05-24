bradley.terry.post=function(theta,data)
{
N=dim(data)[1]; M=length(theta)
sigma=exp(theta[M])
logf=function(k)
{
i=data[k,1]; j=data[k,2]
p=exp(theta[i]-theta[j])/(1+exp(theta[i]-theta[j]))
data[k,3]*log(p)+data[k,4]*log(1-p)
}
sum(sapply(1:N,logf))+sum(dnorm(theta[-M],0,sigma,log=TRUE))
}
