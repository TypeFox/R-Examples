logitnoinflpost=function(beta,y,X)
{
# Log-posterior of beta for the probit model under
# non-informative prior

if (is.matrix(beta)==FALSE) beta=as.matrix(t(beta))
n=dim(beta)[1]
k=dim(beta)[2]
pll=rep(0,n)
for (i in 1:n)
{
lF1=plogis(X%*%beta[i,],log.p=T)
lF2=plogis(-X%*%beta[i,],log.p=T)
pll[i]=sum(y*lF1+(1-y)*lF2)-k/2*log(t(beta[i,])%*%t(X)%*%X%*%beta[i,]/2)+lgamma(k/2)-(k/2)*log(pi)+0.5*log(det(t(X)%*%X))
}
pll
}
