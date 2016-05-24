clme <-
function(dat.obj,theta0, random.names=NULL, no.iter,trace=F,tol.eps=1e-6, method="BFGS")
{
  dat=dat.obj$dat
  lam=dat.obj$lam
  ncoef=length(lam)
  nn=dim(dat)[1]
  ids=unique(dat[,1])
  mm=length(ids)
  num.random=length(random.names)
  cov.names=c(colnames(dat)[1], unique(c(colnames(dat)[2], random.names)))
  cov.index=match(cov.names, colnames(dat), 0L)
  rand.dat=dat[, cov.index]
  tdat.tmp=NULL
  tdat.tmp$random.names=random.names
  
  if(num.random==2){
    random.index=match(random.names, colnames(dat), 0L)
    ss0=theta0[1]
    ss1=theta0[2]
    rho=theta0[3]
    ss=theta0[4]
    ss.vec0=c(ss0,ss1,rho,ss)
    Sigma=matrix(c(ss0^2,rho*ss0*ss1,rho*ss0*ss1,ss1^2),2,2)
    names.error=c("sigma0","sigma1","rho","sigma")
  }else if(num.random==1){
    random.index=match(random.names, colnames(dat), 0L)
    ss.r=theta0[1]
    ss=theta0[2]
    Sigma=ss.r^2
    ss.vec0=c(ss.r, ss)
    names.error=c("sigma.r", "sigma")
  }else if(num.random==0){
    ss=theta0
    ss.vec0=ss
    names.error="sigma"
  }
  
  
  
  beta.vec0=rep(0,ncoef)
  beta.vec=beta.vec0
  
  
  
  ys=rep(0,nn)
  xmats=matrix(0,nrow=nn,ncol=ncoef)
  mat.inv=matrix(0,nrow=nn,ncol=nn)
  
  max.dtheta=1
  i=1

  while(i<=no.iter & max.dtheta>=tol.eps)
  {
    for(j in 1:mm)
    {
      id.x=(dat[,1]==ids[j])
      tmp.y=as.vector(dat[id.x, 3])
      tmp.x=as.matrix(dat[id.x,4:(3+ncoef)])
      if(num.random!=0){
        Z=dat[id.x, random.index]
        Z=as.matrix(Z)
        ww=dim(Z)[1]
        SS=Z%*%Sigma%*%t(Z)+ss^2*diag(ww)
        SS.sqrt=sqrt.mat(SS)
        
        
        SS.inv.sqrt=solve(SS.sqrt)
        mat.inv[id.x,id.x]=SS.inv.sqrt
        ys[id.x]=SS.inv.sqrt%*%tmp.y
        xmats[id.x,]=SS.inv.sqrt%*%tmp.x
      }else{
        ys[id.x]=tmp.y
        xmats[id.x,]=tmp.x
      }
    }
    
    aa1=as.matrix(xmats[,lam==0])
    aa2=as.matrix(xmats[,lam==1])
    aa=cls(y=ys,X=aa2-Px(aa1)%*%aa2)
    
    beta2=aa$betahat
    beta1=solve(t(aa1)%*%aa1)%*%t(aa1)%*%(ys-aa2%*%beta2)
    
    beta.vec[lam==0]=beta1
    beta.vec[lam==1]=beta2
    
    #print("OK")
    
    ww=ys-(aa$yhat+Px(aa1)%*%ys)
    
    
    
    if(num.random!=0){
      ww1=ww
      for(j in 1:mm)
      {
        id.x=(dat[,1]==ids[j])
        SS.inv.sqrt=mat.inv[id.x,id.x]
        SS.sqrt=solve(SS.inv.sqrt)
        ww1[id.x]=SS.sqrt%*%ww[id.x]
      }
      
      yhat=dat[,3]-ww1
      #plot(dat[,"TIME"],ww1,type="p")
      
      #hist(ww1)
      
      
      tdat=cbind(rand.dat, ww=ww1)
      if(num.random==2){
        theta.start=c(log(ss0), log(ss1), log((1+rho)/(1-rho)), log(ss))
      }else if (num.random==1){
        theta.start=c(log(ss.r), log(ss))
      }
      
      
      tdat.tmp$dat=tdat

      
      
      tfit.mle=try(lifetime.mle(dat=tdat.tmp, minusloglik=minus.loglik.lme, starts=theta.start, method),silent=T)
      mle.flag=(attr(tfit.mle,"class")=="try-error")

      

      
      if(length(mle.flag)!=0)
      {
        tfit.mle=try(lifetime.mle(dat=tdat.tmp, minusloglik=minus.loglik.lme, starts=theta.start, method = "Nelder-Mead"),silent=T)
      }

     # print(theta.start)
      mle.flag=(attr(tfit.mle,"class")=="try-error")
      
      if(length(mle.flag)==0)
      {
        tfit=mle.obj.to.fit.obj(obj=tfit.mle)

        
        ss.vec=tfit$ss.vec
        if(num.random==2){
          ss0=ss.vec[1]
          ss1=ss.vec[2]
          rho=ss.vec[3]
          ss=ss.vec[4]
          Sigma=matrix(c(ss0^2,rho*ss0*ss1,rho*ss0*ss1,ss1^2),2,2)
        }else if(num.random==1){
          ss.r=ss.vec[1]
          ss=ss.vec[2]
          Sigma=ss.r^2
        }
        
        max.dtheta=max(abs(c(beta.vec,ss.vec)-c(beta.vec0,ss.vec0)))
        if(trace)
        {
          cat(" criterion=",max.dtheta,"\n")
        }
        
        beta.vec0=beta.vec
        ss.vec0=ss.vec
        
        
        fitted=yhat+tfit$fitted
        i=i+1
      }else if(length(mle.flag)!=0){
        res=list(conv=F)
        return(res)
      }
      
      
    }else if(num.random==0){
      yhat=fitted=ys-ww
      ss.vec=ss=sd(ww)
      error=dat[,3]-fitted
      
      std.error=error/ss
      gen.std.error=ww
      
      
      coef=c(beta.vec,ss.vec)
      names(coef)=c(colnames(dat)[4:(ncoef+3)],names.error)
      loglik=-0.5*nn*log(2*pi)-nn*log(ss)-sum(ww^2)/(2*ss^2)
      tfit=ran.eff=NULL
      res=list(dat=dat.obj,coef=coef,fitted=fitted,yhat=yhat,tfit=tfit,
               error=error,std.error=std.error,loglik=loglik,conv=T,
               gen.std.error=gen.std.error,ran.eff=ran.eff)
      return(res)
    }
        
  }

  error=dat[,3]-fitted
  
  std.error=error/ss
  gen.std.error=ww

  
  coef=c(beta.vec,ss.vec)
  names(coef)=c(colnames(dat)[4:(ncoef+3)],names.error)
  
  #############################################################################

  loglik=-tfit.mle$min
  ran.eff=tfit$coef$random[[1]]

  

  
  res=list(dat=dat.obj,coef=coef,fitted=fitted,yhat=yhat,tfit=tfit,
           error=error,std.error=std.error,loglik=loglik,conv=T,
           gen.std.error=gen.std.error,ran.eff=ran.eff)
  return(res)
}
