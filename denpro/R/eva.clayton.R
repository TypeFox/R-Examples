eva.clayton<-function(x,t,marginal="unif",sig=c(1,1),df=1)
# t>0
{
u<-x[1]
v<-x[2]
marg1<-1
marg2<-1

if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[1])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}
if (marginal=="student"){
   u<-pt(x[1],df)
   v<-pt(x[2],df)
   marg1<-dt(x[1],df)
   marg2<-dt(x[2],df)
}

val<-(1+t)*(u*v)^(-1-t)*(u^(-t)+v^(-t)-1)^(-2-1/t)*marg1*marg2

#if (val<0) val<-0

return(val)
}

