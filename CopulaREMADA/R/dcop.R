
# bivariate Frank copula density
# cpar = copula parameter >0 or <0 (limit case of 0 is independence)
dfrk=function(u,v,cpar)
{ t1=1.-exp(-cpar);
  tem1=exp(-cpar*u); tem2=exp(-cpar*v);
  pdf=cpar*tem1*tem2*t1;
  tem=t1-(1.-tem1)*(1.-tem2);
  pdf=pdf/(tem*tem);
  pdf
}

# Clayton copula density
# cpar = copula parameter >0 
dcln=function(u,v,cpar)
{ u[u==1]<-0.9999999999
  u[u==0]<-0.0000000001 
  v[v==1]<-0.9999999999
  v[v==0]<-0.0000000001 
  tem1=u^(-cpar); tem2=v^(-cpar);
  pdf=(tem1+tem2-1)^(-1/cpar-2)*(1+cpar)*tem1*tem2/(u*v)
  pdf
}

dcln90=function(u,v,cpar)
{ u[u==1]<-0.9999999999
  u[u==0]<-0.0000000001 
  v[v==1]<-0.9999999999
  v[v==0]<-0.0000000001 
  cpar=-cpar
  u=1-u
  dcln(u,v,cpar)
}

dcln180=function(u,v,cpar)
{ u[u==1]<-0.9999999999
  u[u==0]<-0.0000000001 
  v[v==1]<-0.9999999999
  v[v==0]<-0.0000000001 
  u=1-u
  v=1-v
  dcln(u,v,cpar)
}

dcln270=function(u,v,cpar)
{ u[u==1]<-0.9999999999
  u[u==0]<-0.0000000001 
  v[v==1]<-0.9999999999
  v[v==0]<-0.0000000001 
  cpar=-cpar
  v=1-v
  dcln(u,v,cpar)
}




# bivariate normal copula density
# cpar = copula parameter with -1<cpar<1
dbvn=function(u,v,cpar)
{ u[u==1]<-0.9999999999
  u[u==0]<-0.0000000001 
  v[v==1]<-0.9999999999
  v[v==0]<-0.0000000001 
  x1=qnorm(u); x2=qnorm(v)
  qf=x1^2+x2^2-2*cpar*x1*x2
  qf=qf/(1-cpar^2)
  con=sqrt(1-cpar^2)*(2*pi)
  pdf=exp(-.5*qf)/con
  pdf=pdf/(dnorm(x1)*dnorm(x2))
  pdf
}


