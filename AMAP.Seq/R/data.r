RNASeq.Data=function(counts,size=NULL,group,model="nbinom",dispersion=NULL){
  counts=data.matrix(counts)
  group=as.numeric(factor(group))
  nG=nrow(counts)
  if(is.null(size)) size=Norm.GMedian(counts)
  if(is.vector(size)) size=outer(rep(1,nG),size)
  counts=counts[,order(group)]
  size=size[,order(group)]
  group=group[order(group)]
  n1=rowSums(matrix(counts[,group==1],nrow=nG))
  n2=rowSums(matrix(counts[,group==2],nrow=nG))
  s1=rowSums(matrix(size[,group==1],nrow=nG))
  s2=rowSums(matrix(size[,group==2],nrow=nG))

  if(model=="poisson") dispersion=rep(1e-6,nG)
  if(is.null(dispersion))
   if(model=="nbinom"){
    disp.est=dispersion.nb.QL(counts,size,group)
    dispersion=disp.est$dispersion
    lam=sqrt(disp.est$mu1*disp.est$mu2)
    del=log(disp.est$mu1)-log(disp.est$mu2)
   }
  return(list(counts=counts,size=size,group=group,model=model,dispersion=dispersion))
}

Norm.GMedian=function(counts){
  counts=data.matrix(counts)
  counts[counts==0]=.1
  m=apply(counts,1,prod)^(1/ncol(counts))
  s=apply(counts/m,2,median)
  return(s)
}
################
RNASeq.Data.Gene=function(counts,size=NULL,group,model="nbinom",dispersion=NULL,shrink.dispersion=FALSE){
  counts=data.matrix(counts)
  group=as.numeric(factor(group))
  nG=nrow(counts)
  if(is.null(size)) size=apply(counts,2,quantile,.75)
  if(is.vector(size)) size=outer(rep(1,nG),size)
  counts=counts[,order(group)]
  size=size[,order(group)]
  group=group[order(group)]
  n1=rowSums(matrix(counts[,group==1],nrow=nG))
  n2=rowSums(matrix(counts[,group==2],nrow=nG))
  s1=rowSums(matrix(size[,group==1],nrow=nG))
  s2=rowSums(matrix(size[,group==2],nrow=nG))

  if(model=="poisson") dispersion=rep(1e-6,nG)
  if(is.null(dispersion))
   if(model=="nbinom"){
    disp.est=dispersion.nb.QL(counts,size,group)
    dispersion=disp.est$dispersion
    lam=sqrt(disp.est$mu1*disp.est$mu2)
    del=log(disp.est$mu1)-log(disp.est$mu2)
    if(shrink.dispersion)
          dispersion=dispersion.nb.shrink(counts,size,group,lam,del,dispersion)
   }
  return(list(count=counts,size=size,treat=group,model=model,dispersion=dispersion))
}
#
est.v.nb.QL=function(n,mu,v0=NULL){
    n=n+1e-10
    mu=mu+1e-10
    nJ=ncol(n)
    nG=nrow(n)
    nU=100
    if(is.null(v0)) v0=rep(.1,nG)
    u=runif(nU,-2,2)
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

est.mu.nb.mle=function(n,c,v,mu0=NULL){
 if(is.vector(n)) return(n/c)
 if(is.null(mu0)) mu0=rowSums(n)/rowSums(c)
 nJ=ncol(n)
 nG=nrow(n)
 nU=400
 u=runif(nU,-2,2)
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


dispersion.nb.QL=function(counts,size,group,m=NULL,d=NULL,v=NULL,max.iter=5){
   n=data.matrix(counts)+1e-10
   if(is.null(size)) size=apply(counts,2,quantile,.75)
   if(is.vector(size)) size=outer(rep(1,nrow(counts)),size)
   c=data.matrix(size)

   v=rep(.1,nrow(n))
   m1=m2=NULL
   for(i in 1:max.iter){
     m1=est.mu.nb.mle(n[,group==1],c[,group==1],v,mu0=m1)
     m2=est.mu.nb.mle(n[,group==2],c[,group==2],v,mu0=m2)
     mu=cbind(m1,m2)[,group]*c
     v=est.v.nb.QL(n,mu=mu,v0=v)
    }
    return(list(mu1=m1,mu2=m2,dispersion=v))
 }
dispersion.nb.shrink=function(counts,size,group,lam,del,disp){
     v=del
     for(g in 1:nrow(counts)){
       lglk=0
       for(i in 1:length(group)){
        mu=size[g,i]*lam[g]/exp((-1)^group[i]*del[g]/2)
        lglk=lglk+lgamma(counts[g,i]+1/disp)-lgamma(counts[g,i]+1)-lgamma(1/disp)-log(1+mu*disp)/disp-log(1+1/mu/disp)*counts[g,i]
      }
     vlglk=log(disp)+lglk
     v[g]=sum(exp(vlglk-max(vlglk)))/sum(exp(lglk-max(lglk)))*exp(max(vlglk)-max(lglk))
     }
    return(v)
 }

