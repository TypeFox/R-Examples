eva.student<-function(x,t=rep(4,length(x)),
marginal="unif",sig=c(1,1),r=0,df=1)
# t>2 
#  sig is std of marginals
{
if (marginal=="unif"){
   u<-x[1]/sig[1]
   v<-x[2]/sig[2]
   marg1<-1/sig[1]
   marg2<-1/sig[2]
}
if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[1])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}
if (marginal=="student"){
   u<-pt(x[1]/sig[1],df=t[1])
   v<-pt(x[2]/sig[2],df=t[2])
   marg1<-dt(x[1]/sig[1],df=t[1])/sig[1]
   marg2<-dt(x[2]/sig[2],df=t[2])/sig[2]
}

d<-2
x1<-qt(u,df=df)
x2<-qt(v,df=df)
produ<-dt(x1,df=df)*dt(x2,df=df)

nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
vakio<-gamma((df+d)/2)/((df*pi)^(d/2)*gamma(df/2))
ga<-vakio*(1-r^2)^(-1/2)*(1+nelio/df)^(-(df+d)/2)
#ga<-(1-r^2)^(1/2)*(1+(x1^2+x2^2-2*r*x1*x2)/(t*(1-r^2)))^(-(t+d)/2)

val<-ga/produ*marg1*marg2

return(val)
}

