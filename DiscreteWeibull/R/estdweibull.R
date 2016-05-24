estdweibull <-
function(x, method="ML", zero=FALSE, eps=0.0001, nmax=1000)
{
par<-numeric(2)
n<-length(x)
m1<-mean(x)
m2<-mean(x^2)
beta0<-1
q0<-ifelse(zero,m1/(m1+1),(m1-1)/m1)
# method of moments
if(method=="M")
{
if(sum(x<=as.numeric(!zero)+1)==n)
{
message("Method of moments not applicable on this sample!")
par<-c(NA,NA)
}
else
{
par<-solnp(pars<-c(q0,beta0),fun=lossdw,x=x,zero=zero,eps=eps,nmax=nmax,LB=c(0,0),UB=c(1,100))$pars
}
}
# maximum likelihood method
else if (method=="ML")
{
if(sum(x<=as.numeric(!zero)+1)==n)
{
message("Method of maximum likelihood not applicable on this sample!")
par<-c(NA,NA)
}
else
{
par<-nlm(f=loglikedw,x=x,zero=zero,p=c(q0,beta0))$estimate
}
}
# method of proportion
else if (method=="P")
{
y<-sum(x==as.numeric(!zero))
if(y==0)
{
message("Method of proportion not applicable for estimating q!")
par<-c(NA,NA)
}
else
{
par[1]<-1-y/n
z<-sum(x==(as.numeric(!zero)+1))
if(z+y==round(n) | z/n==0)
{
message("Method of proportion not applicable for estimating beta!")
par[2]<-NA
}
else
par[2]<-log(log(par[1]-z/n)/log(par[1]))/log(2)
}
}
par
}

