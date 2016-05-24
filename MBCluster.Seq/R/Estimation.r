cl.mb=function(n,s,t,model,P,C,M,nb.disp=NULL,a=1,iter.max=100){
  C0=C
  if(is.null(P)) P=matrix(1/nrow(C),nrow(n),nrow(C))
  if(is.null(M)){
    M=matrix(0,nrow(n),nrow(C))
    for(k in 1:nrow(C)){
       M[,k]=log(rowSums(n)/rowSums(exp(sweep(s,2,C[k,t],"+"))))
      }
    }
  for(i in 1:iter.max){
    P=cl.mb.est.P(n,s,t,model,P,C,M,nb.disp,a)
    MC=cl.mb.est.MC(n,s,t,model,P,C,M,nb.disp)
    M=MC$M
    C=MC$C   
    dif=apply(abs(C-C0),1,max)/(apply(C,1,max)-apply(C,1,min)+1e-4)
      if(max(dif)<1e-2) break
    C0=C
  }
  print(paste("--------",i,"iterations for a=",a))
  k1=apply(P,1,which.max)
  return(list(P=P,C=C,M=M))
}

###########################################

cl.mb.est.MC=function(n,s,t,model,P,C,M,nb.disp,method="EM"){
  if(method=="SA"){
    P=P==apply(P,1,max)
    P=P+1e-10
  }
  for(k in 1:nrow(C)){
    mc=cl.mb.est.mc(n,s,t,model,p=P[,k],c=C[k,],m=M[,k],nb.disp)
    M[,k]=mc$m
    C[k,]=mc$c
  }
  return(list(C=C,M=M))
}


cl.mb.est.mc=function(n,s,t,model,p,c,m,nb.disp){
  if(model=="poisson"){
     mc=cl.ps.est.mc(n,s,t,p,c)
     m=mc$m
     c=mc$c
   }
  if(model=="nbinom"){
     mc=cl.nb.est.mc(n,s,t,p,c,m,v=nb.disp)
     m=mc$m
     c=mc$c
   }
  return(list(c=c,m=m))
}



##################################

cl.mb.est.P=function(n,s,t,model,P=NULL,C,M,nb.disp,a=1,method="DA"){
 nK=nrow(C)
 lglk=matrix(0,nrow=nrow(n),ncol=nK)
 for(k in 1:nK){
  r=exp(sweep(s+M[,k],2,C[k,t],"+"))
  if(model=="nbinom") 
    lglk[,k]=lglk.nb(n,r,v=nb.disp) 
  if(model=="poisson")  
    lglk[,k]=lglk.ps(n,r) 
 }
 p=rep(1/nK,nK)
 if(!is.null(P)) p=colMeans(P)
 p=matrix(rep(p,each=nrow(n)),nrow=nrow(n))
 f=exp(lglk-maxRow(lglk))
 if(method=="EM") P=p*f/rowSums(p*f)
 if(method=="DA") P=p*f^a/rowSums(p*f^a)
 if(method=="SA") P=p^a*f^a/rowSums(p^a*f^a)
 return(P)
}




cl.nb.est.m=function(n,s,t,c,m=NULL,v){
 n=matrix(n,ncol=length(t))
 s=matrix(s,ncol=length(t))
 v=matrix(v,ncol=length(t))
 r=exp(sweep(s,2,c[t],"+"))
 if(is.null(m)) m=log(rowSums(n)/rowSums(r))
 id=rep(TRUE,nrow(n))
 for(i in 1:20){
  dfm=rowSums(matrix((n[id,]-r[id,]*exp(m[id]))/(1+r[id,]*exp(m[id])*v[id,]),nrow=sum(id)))
  ddfm=rowSums(matrix(-r[id,]*exp(m[id])*(n[id,]*v[id,]+1)/(1+r[id,]*exp(m[id])*v[id,])^2,nrow=sum(id)))
  dm=-dfm/ddfm
  dm[abs(dm)>.5]=sign(dm[abs(dm)>.5])*.5
  m[id]=m[id]+dm
  id[id][abs(dm)<1e-3]=FALSE
  if(sum(id)==0) break
 }
 return(m)
}


cl.nb.est.c=function(n,s,t,p=NULL,c=NULL,m,v){
 n=matrix(n,ncol=length(t))
 s=matrix(s,ncol=length(t))
 if(is.null(p)) p=rep(1,nrow(n))
 if(is.null(c)) c=cl.ps.est.mc(n,s,t,p)$c

 for(i in 1:20){
   r=exp(sweep(sweep(s,1,m,"+"),2,c[t],"+"))
   dfc=sumRow((n-r)*p/(1+r*v),by=t)
   if(nrow(n)>1) dfc=rowSums(t(dfc))
   ddfc=sumRow(-r*(n*v+1)*p/(1+r*v)^2,by=t)
   if(nrow(n)>1) ddfc=rowSums(t(ddfc))
   dc=-dfc/ddfc
   dc[abs(dc)>.5]=sign(dc[abs(dc)>.5])*.5
   if(mean(abs(dc)/(max(c)-min(c)+1e-4))<1e-4) break
   c=c+dc
 }
# print(paste(i,"iterations to estimate c"))
 return(c)
}



cl.nb.est.mc=function(n,s,t,p=NULL,c=NULL,m=NULL,v){
  n=matrix(n,ncol=length(t))
  s=matrix(s,ncol=length(t))
  if(is.null(p))  p=rep(1,nrow(n))
  if(is.null(c)) c=rep(0,length(unique(t)))
  if(is.null(m)){
     r=exp(sweep(s,2,c[t],"+"))
     m=log(rowSums(n)/rowSums(r))
   }
  c0=c
  for(i in 1:100){
    m=cl.nb.est.m(n,s,t,c,m,v)
    c=cl.nb.est.c(n,s,t,p,c,m,v)
    m=m+mean(c)
    c=c-mean(c)
    if(abs(mean(c)-mean(c0))<(max(c)-min(c))*1e-3) break
    c0=c
  }
 return(list(m=m,v=v,c=c))
}


cl.ps.est.mc=function(n,s,t,p=NULL,c0=NULL){
  if(is.null(c0)) c0=rep(0,length(unique(t)))
  n=matrix(n,ncol=length(t))
  s=matrix(s,ncol=length(t))
  if(nrow(n)==1){
    mc=log(sumRow(n,by=t)/sumRow(exp(s),by=t))
    return(list(m=mean(mc),c=mc-mean(mc)))
  }

  if(is.null(p))  p=rep(1,nrow(n))
  n[n<=0]=1e-10
  c=c0

  p=p+(max(p)-min(p))*1e-10
  ps=sumRow(exp(s),by=t)*p
  pn=sumRow(n,by=t)*p
  
  psc=sweep(ps,2,exp(c),"*")
  m=rowSums(pn)/rowSums(psc)
  for( nStep in 1:100){
   psm=ps*m
   c=log(colSums(pn)/colSums(psm))
   c=c-mean(c)
   if(max(abs(c-c0))/(max(c0)-min(c0)+1e-6)<1e-4) break
   c0=c
   psc=sweep(ps,2,exp(c),"*")
   m=rowSums(pn)/rowSums(psc)
  }
  m=log(rowSums(n)/rowSums(exp(sweep(s,2,c[t],"+"))))
  return(list(c=c,m=m))
}




################ estimate dispersion with QL

est.nb.v.QL.one=function(n,mu,v0=NULL){
    n=n+1e-10
    mu=mu+1e-10
    nJ=ncol(n)
    nG=nrow(n)
    nU=400
    if(is.null(v0)) v0=rep(.1,nG)
    u=seq(-2,2,length.out=400)
    d0=rep(Inf,nG)
    id=rep(1,nG)
    for(i in 1:nU){
      v=v0*exp(u[i])
      Qlglk=n*log(n/mu)-(n+1/v)*log((n+1/v)/(mu+1/v))
      d=abs(2*rowSums(Qlglk)-(nJ-1))
      id[d<d0]=i
      d0[d<d0]=d[d<d0]
     }
     v=v0*exp(u[id])
    return(v)
 }

est.nb.mu.mle.one=function(n,c,v,mu0=NULL){
 if(is.vector(n)) return(n/c)
 if(is.null(mu0)) mu0=rowSums(n)/rowSums(c)
 nJ=ncol(n)
 nG=nrow(n)
 nU=400
 u=seq(-2,2,length.out=400)
 lglk0=rep(-Inf,nG)
 id=rep(1,nG)
 for(i in 1:nU){
  m=c*mu0*exp(u[i])
  lglk=rowSums(lgamma(n+1/v)-lgamma(n+1)-lgamma(1/v)-log(1+m*v)/v-log(1+1/m/v)*n)
  id[lglk>lglk0]=i
  lglk0[lglk>lglk0]=lglk[lglk>lglk0]
 }
  mu=mu0*exp(u[id])
  return(mu)
}
est.nb.v.QL=function(Count,Normalizer,Treatment,v=NULL,max.iter=5){
   Count=data.matrix(Count)+1e-10
   Size=exp(Normalizer)
   v=rep(.1,nrow(Count))

   m=sumRow(Count,by=Treatment)/sumRow(Size,by=Treatment)
   for(i in 1:max.iter){
     for(j in 1:ncol(m))
     m[,j]=est.nb.mu.mle.one(Count[,Treatment==j],Size[,Treatment==j],v,mu0=m[,j])
     mu=m[,Treatment]*Size
     v=est.nb.v.QL.one(Count,mu=mu,v0=v)
    }
    return(v)
 }


