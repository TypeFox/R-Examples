estdlaplace2 <-
function(x,method="M",err=0.001,parml=c(exp(-1),exp(-1)))
{
if(method=="M")
{
par<-optim(par=c(exp(-1),exp(-1)),fn=loss,x=x,method="L-BFGS-B",lower=c(err,err),upper=c(1-err,1-err))[[1]]
}
else if(method=="ML")
{
if(sum(x>=0)==length(x) | sum(x<0)==length(x))
{
cat("ML fails!")
par<-c(NA,NA)
}
else
{
par<-optim(par=parml, fn=dlaplacelike2,x=x)[[1]]
}
}
else if(method=="P")
{
r0 <- mean(x==0)
rplus <- mean(x>=0)
if(rplus==r0 | r0==0)
{par<-c(NA,NA)}
else{
par1 <- 1-r0/rplus
par2 <- par1^(rplus/(1-rplus))
par <- c(par1,par2)
}
}
else if(method=="MM")
{
par1 <- mean(x[x>=0])/(1+mean(x[x>=0]))
if(sum(x<0)>0)
{
par2 <- 1-1/abs(mean(x[x<0]))
}
else
{
par2<-NA
}
par <- c(par1,par2)
}
par
}
