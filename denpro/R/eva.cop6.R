eva.cop6<-function(x,t,marginal="unif",sig=c(1,1))
# t in [1,\infty)
{
u<-x[1]
v<-x[2]
marg1<-1
marg2<-1

if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[2])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}

val<-t*(1-u)^(t-1)*(1-v)^(t-1)*
((1-u)^t+(1-v)^t-(1-u)^t*(1-v)^t)^(1/t-2)*
(-(1/t-1)*(1-(1-u)^t)*(1-(1-v)^t)+
 (1-u)^t+(1-v)^t-(1-u)^t*(1-v)^t)*marg1*marg2

#if (val<0) val<-0

return(val)
}

