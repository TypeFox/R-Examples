loglinll=function(beta,y,X)
{

# Log-likelihood of the log-linear model (log-posterior of the log-linear
# model under flat prior)

if (is.matrix(beta)==FALSE) beta=as.matrix(t(beta))
n=dim(beta)[1]
pll=rep(0,n)
for (i in 1:n)
{
lF=exp(X%*%beta[i,])
pll[i]=sum(dpois(y,lF,log=T))
}
pll
}
