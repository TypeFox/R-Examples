tailfunc<-function(R,d,type,gnum=1000,sig=1,nu=1)
{
volball<-function(r,d){ return(r^d*pi^(d/2)/gamma(d/2+1)) }
volsphere<-function(d){ return(2*pi^(d/2)/gamma(d/2)) }

if (type=="bartlett"){
   norma<-d*(d+2)/(2*volsphere(d))
   funni<-function(t,d=d,nu=nu){ return( t^(d-1)*(1-t^2) ) }
   levfun<-function(t,d=d,sig=sig,nu=nu){ return( 1-(t/sig)^2 ) }
}
if (type=="gauss"){
   norma<-(2*pi)^(-d/2)
   funni<-function(t,d=d,nu=nu){ return( t^(d-1)*exp(-t^2/2) ) }
   levfun<-function(t,d=d,sig=sig,nu=nu){ return( exp(-(t/sig)^2/2) ) }
}
if (type=="student"){
   norma<-gamma((nu+d)/2)/((pi*nu)^(d/2)*gamma(nu/2))
   funni<-function(t,d=d,nu=nu){ return( t^(d-1)*(1+t^2/nu)^(-(d+nu)/2) ) }
   levfun<-function(t,d=d,sig=sig,nu=nu){ return( (1+(t/sig)^2/nu)^(-(d+nu)/2) ) }
}

# probability calc (numerical integral)
# y[r] = int_0^(r/sig) funni(t) dt
stepy<-R/sig/gnum
radiy<-seq(stepy,R/sig,stepy)
y<-matrix(0,length(radiy),1)
y[1]<-stepy*funni(radiy[1],d=d,nu=nu)
for (i in 2:length(y)){
    y[i]<-y[i-1]+stepy*funni(radiy[i],d=d,nu=nu)
}

# level calc
step<-R/gnum
radi<-seq(step,R,step)
level<-matrix(0,length(radi),1)
for (i in 1:length(level)){
    level[i]<-levfun(radi[i],d=d,sig=sig,nu=nu)
}

intrad2lev<-0
step<-radi[2]-radi[1]
for (i in 1:length(level)){
   intrad2lev<-intrad2lev+step*level[i]
}
levelnorma<-level/intrad2lev/2

proba<-norma*volsphere(d)*y
volu<-volball(radi,d)
level<-sig^(-d)*norma*level

return(list(radi=radi,proba=proba,volu=volu,level=level,levelnorma=levelnorma))
}









