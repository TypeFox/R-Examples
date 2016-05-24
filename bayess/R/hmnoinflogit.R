hmnoinflogit=function(niter,y,X,scale=length(y))
{
p=dim(X)[2]
mod=summary(glm(y~-1+X,family=binomial(link="logit")))
beta=matrix(0,niter,p)
beta[1,]=as.vector(mod$coeff[,1])
Sigma2=as.matrix(mod$cov.unscaled)
for (i in 2:niter)
{
tildebeta=rmnorm(1,beta[i-1,],scale*Sigma2)
llr=logitnoinflpost(tildebeta,y,X)-logitnoinflpost(beta[i-1,],y,X)
if (runif(1)<=exp(llr)) beta[i,]=tildebeta else beta[i,]=beta[i-1,]
}
beta
}
