loglinnoinflpost=function(beta,y,X)
{

# Log-posterior of beta for the log-linear model under
# non-informative prior

if (is.matrix(beta)==F) beta=as.matrix(t(beta))
n=dim(beta)[1]
k=dim(beta)[2]
pll=rep(0,n)
for (i in 1:n)
{
lF=exp(X%*%beta[i,])
pll[i]=sum(dpois(y,lF,log=T))-(2*k-1)/4*log(t(beta[i,])%*%t(X)%*%X%*%beta[i,]/2)+lgamma((2*k-1)/4)-(k/2)*log(pi)+0.5*log(det(t(X)%*%X))
}
pll
}
