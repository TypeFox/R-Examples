

dst.pairs=function(n,s,t,model,nb.disp,m,C,k,method,nNearest=NULL,pairs=NULL){
 k0=unique(k)
 nK=length(k0)
 if(is.null(pairs)){
     pairs=cbind(rep(k0,each=nK),rep(k0,nK))
     pairs=pairs[pairs[,1]<pairs[,2],]
  }
 if(is.vector(pairs)) pairs=matrix(pairs,1)
 d=rep(0,nrow(pairs))
 for(i in 1:nrow(pairs)){
    n1=n[k==pairs[i,1],]
    n2=n[k==pairs[i,2],]
    s1=s[k==pairs[i,1],]
    s2=s[k==pairs[i,2],]
    v1=nb.disp[k==pairs[i,1],]
    v2=nb.disp[k==pairs[i,2],]
    m1=m[k==paste(pairs[i,1])]
    m2=m[k==paste(pairs[i,2])]
    c1=C[paste(pairs[i,1]),]
    c2=C[paste(pairs[i,2]),]
    if(method=="Ward") d[i]=dst.Ward(n1,n2,s1,s2,m1,m2,c1,c2,v1,v2,t,model)
    if(method=="KL") d[i]=dst.KL(n=n1,s=s1,t,c=c1,c0=c2)/2+dst.KL(n=n2,s=s2,t,c=c2,c0=c1)/2
   } 
  D=cbind(pairs,d)
  D=matrix(D,nrow=nrow(pairs))
  D=D[order(d),]
  D=matrix(D,nrow=nrow(pairs))
  if(nrow(D)>1) D=D[(1:nrow(D))<=nNearest,]
  return(D)
 }
dst.KL=function(n,s,t,c,c0){
  n=matrix(n,ncol=length(t))
  s=matrix(s,ncol=length(t))
  n=n+1e-10
  r=sweep(s,2,c[t],"+")
  r=exp(r)*(rowSums(n)/rowSums(exp(r)))
  r0=sweep(s,2,c0[t],"+")
  r0=exp(r0)*(rowSums(n)/rowSums(exp(r0)))
  d=-r+n*log(r)+r0-n*log(r0)
  d=sum(d)/nrow(n)
}
dst.Ward=function(n1,n2,s1,s2,m1,m2,c1,c2,v1,v2,t,model){
  n1=matrix(n1,ncol=length(t))
  s1=matrix(s1,ncol=length(t))
  n2=matrix(n2,ncol=length(t))
  s2=matrix(s2,ncol=length(t))

  r1=exp(sweep(s1+m1,2,c1[t],"+"))
  r2=exp(sweep(s2+m2,2,c2[t],"+"))
  n=rbind(n1,n2)
  s=rbind(s1,s2)
  v=rbind(v1,v2)
  mc=cl.mb.est.mc(n,s,t,model,p=NULL,c=NULL,m=NULL,nb.disp=v)
  m=mc$m
  c=mc$c
  r=exp(sweep(s+m,2,c[t],"+"))
  if(model=="nbinom") d=sum(lglk.nb(n1,r1,v1))+sum(lglk.nb(n2,r2,v2))-sum(lglk.nb(n,r,v))
  if(model=="poisson") d=sum(lglk.ps(n1,r1))+sum(lglk.ps(n2,r2))-sum(lglk.ps(n,r))
  return(d)
}




