# parameter estimation
estdweibull3<-function(x,method="P",eps=1e-04)
{
n<-length(x)
y<-sum(x==0)
z<-sum(x==1)
# method of proportion
if (method=="P")
{
c<--log(1-y/n)
if(y+z!=n & z!=0 & y!=0)
{
beta<-log(log(1-(y+z)/n)/log(1-y/n)-1)/log(2)
}
else
{
beta<-NA
}
}
# maximum likelihood
else if(method=="ML")
{
c0<-1
beta0<-0
if(sum(x>1)==0)
{
c<-NA
beta<-NA
}
else
{
sol<-optim(p<-c(c0,beta0),fn=loglikedw3,x=x)$par
c<-sol[1]
beta<-sol[2]
}
}
# method of moments
else if(method=="M")
{
if(sum(x>1)==0)
{
c<-NA
beta<-NA
}
else
{
c0<-log((1+mean(x))/mean(x))
beta0<-0
sol<-optim(p<-c(c0,beta0),fn=lossdw3,x=x,eps=eps)$par
c<-sol[1]
beta<-sol[2]
}
}
return(c(hatc=c,hatbeta=beta))
}
