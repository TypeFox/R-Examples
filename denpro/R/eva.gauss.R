eva.gauss<-function(x,t=1,marginal="unif",sig=c(1,1),r=0,tapa1=TRUE)
{
#  sig is std of marginals

if (marginal=="unif"){
   u<-x[1]/sig[1]+1/2
   v<-x[2]/sig[2]+1/2
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
   u<-pt(x[1]/sig[1],df=t)
   v<-pt(x[2]/sig[2],df=t)
   marg1<-dt(x[1]/sig[1],df=t)/sig[1]
   marg2<-dt(x[2]/sig[2],df=t)/sig[2]
}

d<-2
x1<-qnorm(u,sd=1)
x2<-qnorm(v,sd=1)

if (tapa1){
  produ<-dnorm(x1,sd=1)*dnorm(x2,sd=1)
  nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
  vakio<-(2*pi)^(-d/2) 
  g<-vakio*(1-r^2)^(-1/2)*exp(-(1/2)*nelio)
  val<-g/produ*marg1*marg2
}
else{
  # x1,x2 -> copuval
  copuval<-(1-r^2)^(-1/2)*
  exp(-(x1^2+x2^2-2*r*x1*x2)/(2*(1-r^2)))/exp(-(x1^2+x2^2)/2)
  val<-copuval*marg1*marg2
}

return(val)
}

