decom.initial=function(data,nK,p0=NULL,d0=0,nK0=0){
    n1=rowSums(as.matrix(data$counts[,data$group==1]))     
    n2=rowSums(as.matrix(data$counts[,data$group==2]))    
    s1=rowSums(as.matrix(data$size[,data$group==1]))    
    s2=rowSums(as.matrix(data$size[,data$group==2]))    
    keep=n1>0 & n2>0
    n1=n1[keep]
    n2=n2[keep]
    s1=s1[keep]
    s2=s2[keep]
    nG=sum(keep)
    m=sqrt(n1*n2/s1/s2)
    d=log(n1/s1)-log(n2/s2)
    v=rep(1,nG)
    if(nK>1) v=kmeans(cbind(d,d),centers=nK,nstart=10)$cluster
    P=c()  
    for(k in 1:nK) P=cbind(P,v==k)
    MGN=MAX.step(P+1e-10,m,d,p0=NULL,d0=NULL,nK0=0)
    if(nK0==0)     return(MGN)
    pr=MGN$pr
    alpha=MGN$alpha
    beta=MGN$beta
    mu=MGN$mu
    sigma=MGN$sigma

       pr0=alpha0=beta0=mu0=sigma0=rep(0,nK0)
       if(!is.null(p0)) id0=abs(d-d0)<quantile(abs(d-d0),p0)
                   else id0=abs(d-d0)<.1
       p0=mean(id0)
       m0=m[id0]
       v=rep(1,length(m0))
       if(nK0>1) v=sample(nK0,length(m0),replace=TRUE)
       for(k in 1:nK0){
          id=v==k
          pr0[k]=mean(id)
          alpha0[k]=mean(m0[id])^2/var(m0[id])
          beta0[k]=mean(m0[id])/var(m0[id])
          mu0[k]=d0
          sigma0[k]=1e-4
       }
       pr=c(pr*(1-p0),pr0*p0)
       alpha=c(alpha,alpha0)
       beta=c(beta,beta0)
       mu=c(mu,mu0)
       sigma=c(sigma,sigma0)
       MGN=data.frame(pr=pr,alpha=alpha,beta=beta,mu=mu,sigma=sigma)
      return(MGN)
}


################################
EXP.part.poisson=function(n,c,t,a,b,mu,sigma){ 
   n1=sum(n[t==1])
   n2=sum(n[t==2])
   c1=sum(c[t==1])
   c2=sum(c[t==2])
   
   nR=1000  
   r=seq(mu-10*sigma,mu+10*sigma,length.out=nR)
   wd.r = r[-1]-r[-nR]     
   r =  r[-1]/2+r[-nR]/2   
   
   logf.r=-(r-mu)^2/2/sigma^2-log(2*pi*sigma^2)/2+r/2*(n1-n2)-(a+n1+n2)*log(b+c1*exp(r/2)+c2/exp(r/2))+log(wd.r)
   scale=max(logf.r)

   lglk=log(sum(exp(logf.r-scale)))+scale 
   lglk=lglk+a*log(b)-lgamma(a)+lgamma(a+n1+n2)+n1*log(c1)+n2*log(c2)-lgamma(n1+1)-lgamma(n2+1)
   
   lam=sum((a+n1+n2)/(b+c1*exp(r/2)+c2/exp(r/2))*exp(logf.r-scale))/sum(exp(logf.r-scale))
   del=sum(r*exp(logf.r-scale))/sum(exp(logf.r-scale))
   return(c(lglk,lam,del))
}



EXP.part.nbinom=function(n,c,v,t,m,d,log.w=NULL){ 
   if(is.null(log.w)) log.w=-log(length(m))
   lglk=log.w
   for(i in 1:length(t)){
    mu=c[i]*m/exp((-1)^t[i]*d/2)
    lglk=lglk+lgamma(n[i]+1/v)-lgamma(n[i]+1)-lgamma(1/v)-log(1+mu*v)/v-log(1+1/mu/v)*n[i]
    }
   lkhood=exp(lglk-max(lglk))
   lam=sum(m*lkhood)/sum(lkhood)
   del=sum(d*lkhood)/sum(lkhood)
   lglk=log(sum(exp(lglk-max(lglk))))+max(lglk)
   return(c(lglk,lam,del))
}

EXP.step=function(data,MGN,p0,d0,nK0,nMC){
   nG=nrow(data$counts)
   nK=length(MGN$pr)
   P=lglk.gk=lam.gk=del.gk=matrix(0,nG,nK)
   counts=data$counts
   size=data$size
   group=data$group
   dispersion=data$dispersion
   LAM=DEL=matrix(0,nrow=nMC,ncol=nK)
   for(k in 1:nK){
     LAM[,k]=rgamma(nMC,MGN$alpha[k],MGN$beta[k])
     DEL[,k]=rnorm(nMC,MGN$mu[k],MGN$sigma[k])
  }
   for( g in 1:nG)
   for( k in 1:nK){
     if(data$model=="poisson")
       exp.part=EXP.part.poisson(counts[g,],size[g,],group,MGN$alpha[k],MGN$beta[k],MGN$mu[k],MGN$sigma[k])
     if(data$model=="nbinom")
       exp.part=EXP.part.nbinom(counts[g,],size[g,],dispersion[g],group,LAM[,k],DEL[,k],log.w=-log(nMC))
     lglk.gk[g,k]=log(MGN$pr[k])+exp.part[1]
      lam.gk[g,k]=MGN$pr[k]*exp.part[2] 
      del.gk[g,k]=MGN$pr[k]*exp.part[3]
    }
   lam=rowSums(lam.gk)
   del=rowSums(del.gk)
   P=exp(lglk.gk-apply(lglk.gk,1,max))
   P=P/rowSums(P)
   if(nK0>0){
     nK=nK-nK0
     del=(del-d0*sum(MGN$pr[nK+1:nK0]))/sum(MGN$pr[1:nK])
     if(is.null(p0)) p0=1/3
     P0=matrix(P[,nK+1:nK0],ncol=nK0)
     P1=matrix(P[,1:nK],ncol=nK)
     P0=P0/sum(colMeans(P0))*p0
     P1=P1/sum(colMeans(P1))*(1-p0)
     P=cbind(P1,P0)
   }
   #  print(round(cbind(lglk.gk,P),3))
   lglk=sum(P*lglk.gk)

  # print(head(round(cbind(m1=data$n1/data$s1,m2=data$n2/data$s2,lam,del,mu1=lam*exp(del/2),mu2=lam/exp(del/2)),2),10))                 ################ print estimation
   return(list(P=P,lam=lam,del=del,lglk=lglk))
}


MAX.step=function(P,lam,del,p0,d0,nK0){
  nK10=ncol(P)
  pr=colMeans(P)
  alpha=beta=mu=sigma=rep(0,nK10)
  for( k in 1:nK10) {
      z=P[,k]
      ## maximize sum u(-b*lam+a*log.lam+a log.b-lgamma(a))
      x=sum(z*lam)/sum(z)
      y=sum(z*log(lam))/sum(z)
      ##  maximize   -b*x+a*y+a*log(b)-lgamma(a)
      alpha[k]=uniroot(function(a) y+log(a/x)-digamma(a),interval=c(1e-10,1e+10))$root
      beta[k]=alpha[k]/x
    
      mu[k]=sum(z*del)/sum(z)
      sigma[k]=sqrt(sum(z*(del-mu[k])^2)/sum(z))
      if(k>nK10-nK0){
        mu[k]=d0
        sigma[k]=1e-4
      }      
   }
    if(!is.null(p0)) {
       pr[1:(nK10-nK0)]=(1-p0)*pr[1:(nK10-nK0)]/sum(pr[1:(nK10-nK0)])
       pr[1:nK0+nK10-nK0]=p0*pr[1:nK0+nK10-nK0]/sum(pr[1:nK0+nK10-nK0])
     }
  return(data.frame(pr=pr,alpha=alpha,beta=beta,mu=mu,sigma=sigma))
}




MGN.EM=function(data,nK,p0=NULL,d0=0,nK0=0,iter.max=10,print.steps=FALSE,MGN0=NULL,model=NULL,nMC=10000){
   
    if(!is.null(model)) data$model=model
    if(is.null(MGN0))  MGN0=decom.initial(data,nK,p0=p0,d0=d0,nK0=nK0)
    MGN=MGN0
    lglk.old=0

    print("Estimating prior distribution ...")

   for( i in 1:iter.max){
    #   print(paste("Decompostion by EM: step", i))
       if(print.steps) print(as.data.frame(MGN))
       est=EXP.step(data,MGN,p0,d0,nK0,nMC=nMC)
       lglk.new=est$lglk
       lam=est$lam
       del=est$del
       P=est$P
       MGN=MAX.step(P,lam,del,p0,d0,nK0)
   #    if(abs(lglk.old-lglk.new)<nG*1e-3) break
       lglk.old=lglk.new
     }
  #  print(paste("--->> Decompostion ends"))
  #  print(MGN)
    return(list(MGN=MGN,del=del,lam=lam))
}

