# Frank
# cpar = copula parameter: cpar>0 or cpar<0; cpar=0 input will not work
# 1-exp(-cpar) becomes 1 in double precision for cpar>37.4
#CondCinv.frank=function(cpar,p,u)
qcondfrk=function(p,u,cpar)
{ cpar0=exp(-cpar)
  cpar1=1-cpar0
  etem=exp(-cpar*u+log(1./p-1.))   
  tem=1.-cpar1/(etem+1.);
  v=(-log(tem))/cpar
  isinf=is.infinite(v)
  v[isinf]=(-log(cpar0+etem[isinf]))/cpar
  v
} 


#CondCinv.clayton=function(cpar,p,u)
qcondcln=function(p,u,cpar)
{ eta=-cpar/(1+cpar)
  tem=p^eta-1
  tem=tem*(u^(-cpar))+1
  tem^(-1/cpar)
}

qcondbvn<-function(p,u,cpar) 
{ pnorm(sqrt(1-cpar*cpar)*qnorm(p)+cpar*qnorm(u))
}

qcondcln90=function(p,u,cpar)
{ cpar=-cpar
  u=1-u
  qcondcln(p,u,cpar)
}


qcondcln180=function(p,u,cpar)
{ u=1-u
  p=1-p
  1-qcondcln(p,u,cpar)
}


qcondcln270=function(p,u,cpar)
{ cpar=-cpar
  p=1-p
  1-qcondcln(p,u,cpar)
}

