
pdf2grid=function(MGN,mL,mR,dL,dR){
       if(mL>=mR | dL>=dR) return(NULL)
       id=rev(order(MGN$sigma))
       p=MGN$pr[id]
       a=MGN$alpha[id]
       b=MGN$beta[id]
       m=MGN$mu[id]
       s=MGN$sigma[id]
       nK=length(p)

              
       r=c()
       for(k in 1:nK){
         rr=m[k]+s[k]*seq(-5,5,length.out=100)
         r=c(rr,r[r<min(rr) | r>max(rr)])
       }
       r=c(r,seq(dL,dR,length.out=100))     ### grid for del 
       r=sort(unique(r))
       r=r[r>=dL & r<=dR]
       dr=r[-1]-r[-length(r)] 
       r=r[-1]/2+r[-length(r)]/2
   
      u=exp(seq(log(mL),log(mR),.05))       ### grid for lam
      u=exp(u[u>=log(mL) & u<=log(mR)])
      du=u[-1]-u[-length(u)] 
      u=u[-1]/2+u[-length(u)]/2

      nU=length(u)
      nR=length(r)
      u=rep(u,each=nR)
      du=rep(du,each=nR)
      r=rep(r,nU)
      dr=rep(dr,nU)
      
      lglk=matrix(0,length(r),nK)
      for(k in 1:nK)
        lglk[,k]=log(p[k])+dnorm(r,m[k],s[k],log=TRUE)+dgamma(u,a[k],b[k],log=TRUE)+log(dr)+log(du)
      scale=apply(lglk,1,max)
      lglk=log(rowSums(exp(lglk-scale)))+scale
      return(list(m=u,d=r,w=lglk))
}
#######################  MAP TEST   #########################

nb.lglk=function(n,c,v,t,m,d,log.w=0){ 
   lglk=log.w
   for(i in 1:length(t)){
    mu=c[i]*m/exp((-1)^t[i]*d/2)
    lglk=lglk+lgamma(n[i]+1/v)-lgamma(n[i]+1)-lgamma(1/v)-log(1+mu*v)/v-log(1+1/mu/v)*n[i]
    }
   lglk=log(sum(exp(lglk-max(lglk))))+max(lglk)
   return(lglk)
}
test.AMAP=function(data,MGN,del.lim=NULL,FC=NULL,print.steps=FALSE,Integration="MC",nMC=NULL){
  if(!is.null(FC)){
     if(FC<=0){ print("Fold Change to Test Must be a Positive Number");return(0)}
     if(FC<1) FC=1/FC
     if(FC==1) FC=1.01
     if(FC>1) del.lim=c(-log(FC),log(FC))
  }
  counts=data$counts
  size=data$size
  group=data$group
  dispersion=data$dispersion
  nG=nrow(counts)
  s1=s0=rep(0,nG)
  nK=length(MGN$pr)
 # S1=S0=matrix(0,nG,nK)
  
  if(Integration=="MC"){
    if(is.null(nMC)) nMC=50000
    k=sample(nK,nMC,replace=TRUE,prob=MGN$pr)
    lam=rgamma(nMC,MGN$alpha[k],MGN$beta[k])
    del=rnorm(nMC,MGN$mu[k],MGN$sigma[k])
    id0=del>=del.lim[1] & del<=del.lim[2]  
    lam0=lam[id0]
    lam1=lam[!id0]
    del0=del[id0]
    del1=del[!id0]
    lgw0=rep(-log(nMC),sum(id0))
    lgw1=rep(-log(nMC),sum(!id0))
   }
  if(Integration=="grid"){
    s1=rowSums(matrix(data$size[,group==1]))
    s2=rowSums(matrix(data$size[,group==2]))
    n1=rowSums(matrix(data$counts[,group==1]))
    n2=rowSums(matrix(data$counts[,group==2]))
    n1[n1==0]=.1
    n2[n2==0]=.1 
    lam=sqrt(n1/s1*n2/s2)
    del=log(n1/s1)-log(n2/s2)
    lam.min=quantile(lam,.01)/2
    lam.max=quantile(lam,.99)*2
    del.min=quantile(del,.01)/2
    del.max=quantile(del,.99)*2
    grid0=pdf2grid(MGN,mL=lam.min,mR=lam.max,dL=del.lim[1],dR=del.lim[2])
    grid1=pdf2grid(MGN,mL=lam.min,mR=lam.max,dL=del.min,   dR=del.lim[1])
    grid2=pdf2grid(MGN,mL=lam.min,mR=lam.max,dL=del.lim[2],dR=del.max   )
    lam0=grid0$m
    del0=grid0$d
    lgw0=grid0$w
    lam1=c(grid1$m,grid2$m)
    del1=c(grid1$d,grid2$d)
    lgw1=c(grid1$w,grid2$w)
  }
  print(paste(length(lgw0)," | ",length(lgw1)+length(lgw0)," points were drawn"))
  s1=s0=rep(0,nG)
  for( g in 1:nG){
       if(print.steps) if(g %in% round(seq(1,nG,length.out=100)))
       print(paste("--> ",floor(g/nG*100),"% finished to test ",round(del.lim[1],3),"< log-FC <",round(del.lim[2],3),sep=""))
       s1[g]=nb.lglk(counts[g,],size[g,],dispersion[g],group,lam1,del1,lgw1)
       s0[g]=nb.lglk(counts[g,],size[g,],dispersion[g],group,lam0,del0,lgw0)
      }   
   scale=s0*(s0>s1)+s1*(s0<=s1)
   s=log(exp(s0-scale)+exp(s1-scale))+scale
   lglk=s
   s=s0-s
   q=p=exp(s)
   od=order(s)
   q[od]=cumsum(q[od])/(1:nG)
  return(data.frame(stat=s,prob=exp(s),fdr=q))
}

