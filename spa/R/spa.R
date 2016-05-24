"spa" <-
  function(y,x,graph,type=c("soft","hard"),kernel=function(r,lam=0.1){exp(-r/lam)},
           global,control,...){
  cl <- match.call()
  cl[[1]]<-as.name("spa")
  if(is.null(y)|missing(y)){
    stop(paste("Error in y:  must not be NULL or missing"))
  }

  if(!missing(graph)){
    if(is.data.frame(graph))
      graph=as.matrix(graph)
    if(!is.matrix(graph))
      stop(paste("Error in graph:  must be of type  matrix"))
  }

  if(missing(graph))
    stop(paste("Error in graph:  must be supplied"))
  
  if(missing(y)){
    stop("Error in y:  a response (either continuous or binary) must be provided")
  }
  xdat=FALSE
  n=dim(graph)[1]
  m=length(y)
  p=0
  if(!missing(x)){
    x=as.matrix(x)
    if(!is.matrix(x)){
     stop("Error in x: can not convert x to a mtrix")
    }
    if(dim(x)[1]!=n)
      stop("Error in dims x: dim(x)[1]!=dim(graph)[1]")
    xdat=TRUE
  }
  if(m<=n){
    if(m<n){
      y1=rep(NA,n)
      y1[1:m]=y
      y=y1
    } 
    L=c(which(!is.na(y)))
    U=setdiff(1:n,L)
    m=length(L)
  }
  if(m>n){
    stop("Error: |y|> dim(graph)[1]")
  }
  ord=c(L,U)
  lev=NULL
  if(missing(control))
    control=spa.control()
  if(missing(type)){
    type="soft"
    if(is.factor(y)){
      if(nlevels(y)>2){
        ##y=sapply(1:nlevels(y),function(i)as.numeric(y==i))
        stop(paste("Error in y: Currently the response must be binary"))
      }
      y=c(0,1)[as.numeric(y)]
    }
  }
  if( (type=="hard") & xdat){
    if(control$warn)
      warning("x ignored in hard labeled spa")
    xdat=FALSE
  }
  
  if(type=="hard"){
    lev<-paste(c(0,1))
    if(is.factor(y))
      lev<-levels(y)
    if(nlevels(as.factor(y))>2)
      stop(paste("Error in y: Currently the response must be binary"))
    y=as.numeric(as.factor(y))
    y=c(0,1)[y]
  }
  if(type!="soft" & type != "hard"){
    warning("type does not make sense:  default=soft (type should be soft or hard)")
    type="soft"
  }
  resp<-as.numeric(y)[L]
  if(missing(global))
    global=mean(resp)
  if(!control$dissimilar)
    kernel=function(x,lam)x
  if(is.null(control$lval)){
    parm.est=spa.gcv(graph,resp,type,lev,kernel,control,L,U,...)
  }else{
    lval=as.numeric(control$lval)
    if(lval<0)
      lval=-lval
    parm.est=list(cvlam=lval)
  }
  eps=control$adjust
  
  W=kernel(graph,parm.est$cvlam)+eps  
  alpL= sum(as.vector(apply(W[L,L],1,sum)/apply(W[L,],1,sum)))/length(L)
  alpU=sum(as.vector(apply(W[U,U],1,sum)/apply(W[U,],1,sum)))/length(U)

  meas=c(alpL,alpU)

  fit=rep(0,n)
  if(xdat){
    p=dim(x)[2]
    ignore=rep(FALSE,n)
    if(length(U)==1){
      ignore[U]=sum(W[U,L])==0
    }else{
      ignore[U]=apply(W[U,L],1,sum)==0
    }
    U1=setdiff(which(!ignore),L)
    n1=length(L)+length(U1)
    tmp=W[c(L,U1),c(L,U1)]
    tmp=tmp/apply(tmp,1,sum)
    tmp=tmp-t(matrix(apply(tmp,2,sum),n1,n1))/n1
    xstr=ftf.x.fit(resp,x[c(L,U1),],tmp,control=control)
    tdf=xstr$tdf
    fit=rep(global,n)
    fit[c(L,U1)]=xstr$yhat
    xstr=list(coefs=xstr$coefs,fg=xstr$fg)
  }else{
    Sj=W/apply(W,1,sum)
    vals=ftf.fit(resp,Sj[c(L,U),c(L,U)],type,control=control,GCV=TRUE)
    fit[c(L,U)]<-as.vector(vals[[1]])
    tdf<-vals[[2]]
    xstr=NULL
  }

  res=resp-fit[L]
  if(type=="soft"){
    err=sum( (res)^2)
    gcv=err/(1-tdf/m)^2
    sigh=sqrt(err/(n-tdf))
    p1<-length(levels(as.factor(y)))==2
    p2<-control$pce
    adjust=p1&p2
    if(adjust){
      fit[fit<0.001]=0.001
      fit[fit>0.999]=0.999
    } 
    conf=c(sqrt(err),gcv,tdf,sigh)
    model=list(fit=fit,res=res,measure=meas,conf=conf,
      parm.est=parm.est,xstr=xstr,dims=c(n,m,p),dat=c(0,Inf),
      y=y[L],L=L,U=U)
  }else{
    K=length(lev)
    probs=matrix(0,n,K)
    probs[c(L,U),]<-cbind(1-fit,fit)
    cls=lev[apply(probs,1,which.max)]
    conf=table(as.factor(cls[L]),as.factor(lev[resp+1]),
      dnn=c("True value","Final Prediction"))
    model=list(fit=fit,cls=cls,measure=meas,probs=probs,lev=lev,
      res=res,conf=conf,parm.est=parm.est,dims=c(n,m,p),
      dat=c(0,Inf),y=y[L],L=L,U=U)
  }
  obj<-structure(list(call=cl,kernel=kernel,type=type,global=global,
                      control=control,model=model),class="spa")
  obj
}

