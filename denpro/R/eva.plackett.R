eva.plackett<-function(x,t,marginal="unif",sig=c(1,1))
# t>=0, t \neq 1
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

val<-t*(1+(t-1)*(u+v-2*u*v))*((1+(t-1)*(u+v))^2-4*t*(t-1)*u*v)^(-3/2)*marg1*marg2
#if (val<0) val<-0

return(val)
}

