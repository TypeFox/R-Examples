factest=function(object,start=NULL,se=F, control=list()){
 N=object$N
 nit=object$nit
 rt=object$rt
 x=object$score
 con=object$control
 model=object$model
 A=object$A
 W=object$W
 nq=length(A)
 ai=object$par.log[1:nit]
 vi=object$par.log[(nit+1):(2*nit)]
 ter=object$par.log[(2*nit+1):(3*nit)]
 sd2A=object$par.log[3*nit+1]
 sd2V=object$par.log[3*nit+2]
 if(!is.null(start) & length(start)!=(2*N))
   stop("If you provide starting values, this vector should be of size 2 x [number of subjects].\n")
 con=list(method="BFGS",eps=.01,delta=1e-2,trace=0,fnscale=1,parscale=rep(1,2*N),maxit=19999,reltol=sqrt(.Machine$double.eps))
 con[names(control)]=control
 if(!is.logical(se)) stop("the 'se' argument should be of type logical")
 
  ou=matrix(,N,2)
  for(i in 1:N){
    xx=calcEZ(mean(x[i,],na.rm=T), var(rt[i,],na.rm=T), mean(rt[i,],na.rm=T), N=ncol(x))
    ou[i,1]=xx[[1]]
    ou[i,2]=xx[[2]]
  }
  ou[is.nan(ou)]=NA
  ou[ou<0]=.1
  ou[,2]=log(ou[,2])
  if(model==2) ou[,1]=log(ou[,1])
   apstart=scale(ou[,2])*sqrt(exp(sd2A))
   vpstart=scale(ou[,1])*sqrt(exp(sd2V))

 rt[is.na(rt)]=-999
 rt[is.na(x)]=-999
 x[is.na(x)]=-999

   pars=c(apstart,vpstart)
  if(!is.null(start)) pars[!is.na(start)]=start[!is.na(start)]

   resA=optim(pars[1:N],fn=LLscoreA,
      rt=rt,x=x,N=N,nit=nit,ai=ai,vi=vi,ter=ter,sd2A=sd2A,sd2V=sd2V,model=model,A=A,W=W,nq=nq,eps=con$eps,delta=con$delta,
      method="BFGS",control=list(trace=con$trace,fnscale=con$fnscale,maxit=con$maxit,reltol=con$reltol),hessian=se)
   if(resA$convergence!=0) stop("Convergence problems. Optim error code: ",resA$convergence,"\n") else parA=exp(resA$par)

  if(model==1){
   resV=optim(pars[(N+1):(2*N)],fn=LLscoreVd,
    rt=rt,x=x,N=N,nit=nit,ai=ai,vi=vi,ter=ter,sd2A=sd2A,sd2V=sd2V,model=model,A=A,W=W,nq=nq,eps=con$eps,delta=con$delta,
    method="BFGS",control=list(trace=con$trace,fnscale=con$fnscale,maxit=con$maxit,reltol=con$reltol),hessian=F)
  if(resV$convergence!=0) stop("Convergence problems. Optim error code: ",resV$convergence,"\n") else parV=resV$par
  }
  if(model==2){
   resV=optim((pars[(N+1):(2*N)]),fn=LLscoreVq,
   rt=rt,x=x,N=N,nit=nit,ai=ai,vi=vi,ter=ter,sd2A=sd2A,sd2V=sd2V,model=model,A=A,W=W,nq=nq,eps=con$eps,delta=con$delta,
   method="BFGS",control=list(trace=con$trace,fnscale=con$fnscale,maxit=con$maxit,reltol=con$reltol),hessian=F)
  if(resV$convergence!=0) stop("Convergence problems. Optim error code: ",resV$convergence,"\n") else parV=exp(resV$par)
  }   
  ret=matrix(c(parA,parV),N,2)
  dimnames(ret)[[2]]=c("gamma[p]","theta[p]")
  if(se) {
    varA=matrix((diag(solve(resA$hessian))),N,1)
    varV=matrix((diag(solve(resV$hessian))),N,1)
    var=cbind(varA,varV) 
    var[,1]=ret[,1]^2*var[,1]
    if(model==2) var=ret^2*var
    dimnames(var)[[2]]=c("gamma[p].se","theta[p].se")
    ret=cbind(ret[,1],sqrt(var[,2]),ret[,2],sqrt(var[,2]))
  }
  return(ret)
  }
 
  