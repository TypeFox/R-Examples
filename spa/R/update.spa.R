"update.spa" <-
  function(object,ynew,xnew,gnew,type=c("vector","probs","coef","all"),
           reg=c("ridge","hlasso"),trans.update=FALSE,
           dat=list(k=0,l=Inf),verbose=FALSE,FUN=sum,...){
    if(class(object)!="spa"){
      stop("Error: object must be of type spa")
    }
    if(missing(type)){
      type="vector"
    }
    if(missing(reg)){
      reg="ridge"
    }
    if(missing(ynew)){
      stop("Error in ynew:  a response (either continuous or binary) must be provided")
    }
    if(missing(gnew))
      stop(paste("Error in graph:  must be supplied"))
    
    if(!missing(gnew)){
      if(is.data.frame(gnew))
        gnew=as.matrix(gnew)
      if(!is.matrix(gnew))
        stop(paste("Error in graph:  must be of type  matrix"))
    }
    n=dim(gnew)[1]
    if(!missing(xnew)){
      xnew=as.matrix(xnew)
      if(dim(xnew)[1]!=n)
        stop("Error in dims xnew: dim(xnew)[1]!=dim(gnew)[1]")
    }
    xstr=NULL
    xdat=!is.null(object$model$xstr)
    if(missing(xnew)& xdat){
      stop("Error: object xnew must be supplied whenever object was fit with x data")
    }
    if(!xdat){
      if(type=="coef")
        return(NULL)
    }
    y=ynew
    m=length(ynew)
    if(m<=n){
      if(m<n){
        y1=rep(NA,n)
        y1[1:m]=ynew
        ynew=y1
      } 
      L=c(which(!is.na(ynew)))
      U=setdiff(1:n,L)
      m=length(L)
    }
    if(m>n){
      stop("Error: |y|> dim(gnew)[1]")
    }
    if(object$type=="hard"){
      if(is.factor(ynew)){
        lev=c(0,1)
      }
      if(any(levels(ynew[L])!=lev)){
        stop(paste("Error in y: The response levels must correspond to object"))
      }
    }

    ynew=as.numeric(as.factor(ynew))
    ynew=c(0,1)[as.numeric(ynew)]
    
    reg=c(0,1)[match(reg,c("ridge","hlasso"))]
    W=object$kernel(gnew,object$model$parm.est$cvlam)+object$control$adjust
    
    wsub=rep(Inf,n)
    wsub[U]=m-apply(W[U,L],1,FUN)
    ctype=object$type
    prior<-object$global
    control=object$control
    
    pen=rep(0,n)
    k<-dat$k
    Vr=max(wsub[U])+1
    estk=1
    if(k>1){
      Vr=unique(quantile(wsub[U],1:(k-1)/(k-1)))
      estk<-length(Vr)
      Vr[estk]=Inf
    }
    vec<-1/(dat$l+1)* 1:estk/estk
    mx=length(Vr)
    
    resp=rep(0,n)
    resp[L]=ynew[L]
    fit=rep(0,n)
    cfit=rep(0,n)
    prn=rep(0,n)
    if(ctype=="hard"){
      lev<-c(0,1)
    }else{
      reg=0
    }
    tL<-tLU<-L  
    for(r in 1:mx){
      tL=tLU
      wsub[tLU]=Inf
      tU=setdiff(which(wsub<=Vr[r]),tL)
      if(length(tU)==0)
        next
      tLU<-c(tL,tU)
      resp[L]<-ynew[L]
      yval<-resp[tL]
      if(xdat){
        n1=length(tLU)
        prn[tU]=rep(prior,length(tU))
        ignore=rep(TRUE,n1)
        ignore[tL]=FALSE
        if(length(tU)==1){
          ignore[tU]=sum(W[tU,tL])==0
        }else{
          ignore[tU]=apply(W[tU,tL],1,sum)==0
        }
        U1=setdiff(which(!ignore),tL)
        n2=length(tL)+length(U1)
        tmp=W[c(tL,U1),c(tL,U1)]
        tmp=tmp/apply(tmp,1,sum)
        tmp=tmp-t(matrix(apply(tmp,2,sum),n2,n2))/n2
        xstr=ftf.x.fit(yval,xnew[c(tL,U1),],tmp,control=control)
        prn[c(tL,U1)]=xstr$yhat
        xstr=list(coefs=xstr$coefs,fg=xstr$fg,cv=xstr$cv)      
      }else{
        smo<-W[tLU,tLU]
        smo<-smo/apply(smo,1,sum)
        vals=ftf.fit(yval,smo,type=ctype,global=prior,control,GCV=TRUE)
        prn[tLU]=vals[[1]]
        tdf=vals[[2]]
      }
      h=vec[r]*(1-reg)+reg*vec[r]/(  sqrt(prior*(1-prior))*(1-vec[r])+vec[r])
      pen[tU]=h
      nprob=(1-h)*prn[tU] +prior*h
      resp[tU]<-nprob
      fit[tU]<-nprob
      
      if(ctype=="class") {
        resp[tU]<-round(resp[tU])
      }
      if(verbose){
        cat("iter=",r,"\t|tU|=",length(tU),"\t|tL|=",length(tL))
        if(xdat)
          cat("\t||B||_2=",sum( (xstr$coef)^2))
        cat("\n") 
      }
    }
    fit[L]=prn[L]
    if(!trans.update){
      if(ctype=="soft"){
        if(type=="vector")
          return(fit)
        if(type=="coef")
          return(xstr$coef)
        return(list(class=fit,pen=pen))
      }
      ty<-as.factor(fit>0.5)
      ty<-as.factor(lev[ty])
      if(type=="vector"){
        return(ty)
      }
      fit[fit>1]<-1
      fit[fit<0]<-0
      if(type=="probs"|type=="prob"){
        return(cbind(1-fit,fit))
      }
      return(list(class=ty,probs=cbind(1-fit,fit),pen=pen))
    }
    alpL=sum(W[L,L])/sum(W[L,])
    alpU=sum(W[U,U])/sum(W[U,])
    meas=c(alpL,alpU)

    res=ynew[L]-fit[L]
    if(ctype=="soft"){
      err=sum( (res)^2)
      p1<-length(levels(as.factor(ynew)))==2
      p2<-object$control$pce
      adjust=p1&p2
      if(adjust){
        fit[fit<0.001]=0.001
        fit[fit>0.999]=0.999
      } 
      conf=object$model$conf
      conf[1]=sqrt(err)
      model=list(fit=fit,res=res,measure=meas,
        conf=conf,parm.est=object$model$parm.est,xstr=xstr,
        dims=c(n,m,object$model$dims[3]),dat=c(dat[[1]],dat[[2]]),y=ynew[L],L=L,U=U)
    }else{
      K=length(lev)
      probs=matrix(0,n,K)
      probs[c(L,U),]<-cbind(1-fit,fit)
      cls=lev[apply(probs,1,which.max)]
      conf=table(as.factor(cls[L]),as.factor(y[L]),
        dnn=c("True value","Final Prediction"))
      model=list(fit=fit,cls=cls,measure=meas,
        probs=probs,lev=lev,res=res,conf=conf,
        parm.est=object$model$parm.est,
        dims=c(n,m,object$model$dims[3]),dat=c(dat[[1]],dat[[2]]),
        y=ynew[L],L=L,U=U)
    }
    object$model=model
    object
  }

