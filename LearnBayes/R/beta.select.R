beta.select=function(quantile1,quantile2)
{
betaprior1=function(K,x,p)
# suppose one is given a beta(K*m, K*(1-m)) prior 
# where the pth quantile is given by x
# function outputs the prior mean m
{
m.lo=0; m.hi=1; flag=0
while(flag==0)
{
m0=(m.lo+m.hi)/2
p0=pbeta(x,K*m0,K*(1-m0))
if(p0<p) m.hi=m0 else m.lo=m0
if(abs(p0-p)<.0001) flag=1
}
return(m0)
}

p1=quantile1$p; x1=quantile1$x
p2=quantile2$p; x2=quantile2$x

logK=seq(-3,8,length=100); K=exp(logK)
m=sapply(K,betaprior1,x1,p1)

prob2=pbeta(x2,K*m, K*(1-m))
ind=((prob2>0)&(prob2<1))
app=approx(prob2[ind],logK[ind],p2)
K0=exp(app$y)
m0=betaprior1(K0,x1,p1)

return(round(K0*c(m0,(1-m0)),2))
}