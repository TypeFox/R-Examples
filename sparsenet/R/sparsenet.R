sparsenet=function(x,y,weights,exclude,dfmax=nvars+1,pmax=min(dfmax*2,nvars),ngamma=9,nlambda=50,max.gamma=150,min.gamma=1.000001,lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4),lambda=NULL,gamma=NULL,parms=NULL,warm=c("lambda","gamma","both"),thresh=1e-5,maxit=1000000){
    this.call=match.call()
    warm=match.arg(warm)
    istart=switch(warm,
      both=3,
      gamma=2,
      lambda=1
      )
    storage.mode(istart)="integer"
    np=dim(x)
    nobs=as.integer(np[1])
    nvars=as.integer(np[2])
    vnames=colnames(x)
    if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")
    storage.mode(x)="double"
    storage.mode(y)="double"
    if(missing(weights))weights=rep(1,nobs)
    storage.mode(weights)="double"
    if(!missing(exclude)){
        jd=match(exclude,seq(nvars),0)
        if(!all(jd>0))stop("Some excluded variables out of range")
        jd=as.integer(c(length(jd),jd))
    }else {
      jd=as.integer(0)
      exclude=NULL
    }
    max.lambda=lambda0(x,y,weights,exclude)
    max.gamma=as.double(max.gamma)
    ne=as.integer(dfmax)
    nx=as.integer(pmax)
### We handle all the parameters, which we hand in via parms
    flmin=as.double(lambda.min.ratio)
    if(is.null(parms)){
          if(is.null(lambda)){
            if(lambda.min.ratio>=1)stop("lambda.min.ratio should be less than 1")
            lambda=exp(seq(from=log(max.lambda),to=log(max.lambda*lambda.min.ratio),length=nlambda))
          }
          else{
            lambda=rev(sort(lambda))
            nlambda=length(lambda)
          }
          if(is.null(gamma)){
            gamma=exp(seq(from=log(max.gamma),to=log(min.gamma),length=ngamma-1))
            gamma=c(9.9e35,gamma)
          }
          else {
            gamma=rev(sort(gamma))
            ngamma=length(gamma)
          }
        parms=array(0,c(2,ngamma,nlambda),dimnames=list(c("gamma","lambda"),paste("g",seq(ngamma),sep=""),paste("l",seq(nlambda),sep="")))
          parms[1,,]=matrix(rep(gamma,rep(nlambda,ngamma)),ngamma,nlambda,byrow=TRUE)
          parms[2,,]=matrix(rep(lambda,ngamma),ngamma,nlambda,byrow=TRUE)
        }
    else
    {
      dd=dim(parms)
      ngamma=dd[2]
      nlambda=dd[3]
    }
    storage.mode(parms)="double"
    igrid=as.integer(1)
    ngamma=as.integer(ngamma)
    nlambda=as.integer(nlambda)

    thresh=as.double(thresh)
    maxit=as.integer(maxit)
    fit=.Fortran("sparsenet",
       nobs,nvars,x,y,weights,jd,ne,nx,ngamma,nlambda,max.gamma,flmin,parms=parms,igrid,istart,thresh,maxit,
       lmu=integer(1),#actual number of lambdas used
       a0=matrix(double(ngamma*nlambda),ngamma,nlambda),
       ca=array(double(nx*ngamma*nlambda),c(nx,ngamma,nlambda)),
       ia=integer(nx),
       nin=matrix(integer(ngamma*nlambda),ngamma,nlambda),
       rsq=matrix(double(ngamma*nlambda),ngamma,nlambda),
       nlp=integer(1),
       jerr=integer(1),
       PACKAGE="sparsenet")
    lmu=fit$lmu
    coeflist=getcoef_list(fit,nvars,nx,vnames,ngamma)
    rsq=fit$rsq[,seq(lmu)]
    dimnames(rsq)=list(paste("g",seq(ngamma),sep=""),paste("l",seq(lmu),sep=""))
    parms=fit$parms[,,seq(lmu)]
    outlist=list(call=this.call,rsq=rsq,jerr=fit$jerr,coefficients=coeflist,parms=parms,gamma=parms[1,,1],lambda=parms[2,1,],max.lambda=max.lambda)
    class(outlist)="sparsenet"
    outlist
}


