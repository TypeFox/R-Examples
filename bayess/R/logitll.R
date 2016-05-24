logitll=function(beta,y,X)
{
# Log-likelihood of the logit model (log-posterior of the probit
# model under flat prior)

if (is.matrix(beta)==F) beta=as.matrix(t(beta))
n=dim(beta)[1]
pll=rep(0,n)
for (i in 1:n)
{
lF1=plogis(X%*%beta[i,],log.p=T)
lF2=plogis(-X%*%beta[i,],log.p=T)
pll[i]=sum(y*lF1+(1-y)*lF2)
}
pll
}
