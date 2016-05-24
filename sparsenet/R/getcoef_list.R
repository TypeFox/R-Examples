getcoef_list=function(fit,nvars,nx,vnames,ngamma){
    lmu=fit$lmu
    if(lmu<1)stop("an empty model has been returned; something is wrong")
    stepnames=paste("l",seq(lmu)-1,sep="")
    gamnames=paste("g",seq(ngamma),sep="")
    coeflist=as.list(gamnames)
    names(coeflist)=gamnames
    betalist=coeflist
    dd=c(nvars,lmu)
    nins=fit$nin[,seq(lmu)]
    a0s=fit$a0
    parms=fit$parms
    ninmax=max(nins)
    if(ninmax>0){
      cas=array(fit$ca[seq(nx*ngamma*lmu)],c(nx,ngamma,lmu))[seq(ninmax),,,drop=FALSE]
      ja=fit$ia[seq(ninmax)]
      oja=order(ja)
      ja=rep(ja[oja],lmu)
      ia=cumsum(c(1,rep(ninmax,lmu)))
### See getcoef in glmnet package for what is going on here
        for(i in seq(ngamma)){
          ca=cas[,i,]
          df=apply(abs(ca)>0,2,sum)
          beta=new("dgCMatrix",Dim=dd,Dimnames=list(vnames,stepnames),x=as.vector(ca[oja,]),p=as.integer(ia-1),i=as.integer(ja-1))
          betalist[[i]]=list(beta=drop0(beta),df=df)
         }
    }else
      for(i in seq(ngamma)){
        beta = zeromat(nvars,lmu,vnames,stepnames)
        df=rep(0,lmu)
        betalist[[i]]=list(beta=beta,df=df)
      }
  for(i in seq(ngamma)){
    lambda=parms[2,i,seq(lmu)]
    gamma=parms[1,i,1]
    a0=a0s[i,seq(lmu)]
    names(a0)=stepnames
    coeflist[[i]]= c(betalist[[i]],list(a0=a0,dim=dd,lambda=lambda,gamma=gamma))
    }
coeflist
}
