eva.copula<-function(x,type="gauss",marginal="unif",sig=rep(1,length(x)),r=0,
t=rep(4,length(x)),g=1)
{
# sig is std of marginals, r is the correlation coeff, 
# t is deg of freedom

d<-length(x)
marg<-matrix(0,d,1)
u<-matrix(0,d,1)

if (marginal=="unif"){
   for (i in 1:d){
      u[i]<-x[i]/sig[i]  #+1/2
      marg[i]<-1/sig[i]
   }
}
if ((marginal=="normal")||(marginal=="gauss")){
   for (i in 1:d){
      u[i]<-pnorm(x[i]/sig[i])
      marg[i]<-evanor(x[i]/sig[i])/sig[i]
   }
}
if (marginal=="student"){
   for (i in 1:d){
      u[i]<-pt(x[i]/sig[i],df=t[i])
      marg[i]<-dt(x[i]/sig[i],df=t[i])/sig[i]
   }
}
if (type=="gauss"){
   d<-2
   x1<-qnorm(u[1],sd=1)
   x2<-qnorm(u[2],sd=1)

#   produ<-dnorm(x1,sd=1)*dnorm(x2,sd=1)
#   nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
#   vakio<-(2*pi)^(-d/2) 
#   g<-vakio*(1-r^2)^(-1/2)*exp(-(1/2)*nelio)
#   val<-g/produ*marg[1]*marg[2]

  copuval<-(1-r^2)^(-1/2)*
  exp(-(x1^2+x2^2-2*r*x1*x2)/(2*(1-r^2)))/exp(-(x1^2+x2^2)/2)
  val<-copuval*marg[1]*marg[2]

}
if (type=="gumbel"){
  link<-function(y,g){ return ( (-log(y))^g ) }
  linkinv<-function(y,g){ return ( exp(-y^(1/g)) ) }
  der1<-function(y,g){ return ( -g*(-log(y))^(g-1)/y ) }
  der2<-function(y,g){ return ( g*y^(-2)*(-log(y))^(g-2)*(g-1-log(y)) ) }
  
  linky<-link(u,g)
  a<-sum(linky)
  b<-linkinv(a,g)
  der1b<-der1(b,g)
  der2b<-der2(b,g)
  psi<--der2b*der1b^(-3)
  deriy<-der1(u,g) 
  val<-psi*abs(prod(deriy))*prod(marg)
}
if (type=="frank"){
  link<-function(y,g){ return ( -log((exp(-g*y)-1)/(exp(-g)-1)) ) }
  linkinv<-function(y,g){ return ( -log(1+(exp(-g)-1)/exp(y))/g  ) }
  der1<-function(y,g){ return ( g*exp(-g*y)/(exp(-g*y)-1) ) }
  der2<-function(y,g){ return ( g^2*exp(-g*y)/(exp(-g*y)-1)^2 ) }
  
  linky<-link(u,g)
  a<-sum(linky)
  b<-linkinv(a,g)
  der1b<-der1(b,g)
  der2b<-der2(b,g)
  psi<--der2b*der1b^(-3)
  deriy<-der1(u,g) 
  val<-psi*abs(prod(deriy))*prod(marg)
}
if (type=="clayton"){
  link<-function(y,g){ return ( y^(-g)-1 ) }
  linkinv<-function(y,g){ return ( (y+1)^(-1/g) ) }
  der1<-function(y,g){ return ( -g*y^(-g-1) ) }
  der2<-function(y,g){ return ( g*(g+1)*y^(-g-2) ) }
  
  linky<-link(u,g)
  a<-sum(linky)
  b<-linkinv(a,g)
  der1b<-der1(b,g)
  der2b<-der2(b,g)
  psi<--der2b*der1b^(-3)
  deriy<-der1(u,g) 
  val<-psi*abs(prod(deriy))*prod(marg)
}

return(val)
}






