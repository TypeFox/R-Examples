
print.hmmspec <- function(x, ...){
  cat("Hidden Markov Model specification:\n")
  cat(sprintf("J (number of states): \n%i \n", x$J))
  cat("init:\n")
  print(x$init)
  cat ("transition:\n")
  print(x$transition)
  cat("emission:\n")
  print(x$parms.emission)
  return(invisible(x))
}

summary.hmm <- function (object, ...) 
{
    cat("init: \n", round(object$model$init, 3), "\n")

    cat ("\ntransition:\n")
    print(round(object$model$trans, 3)) ## FIXME: This is **Very** fragile
    cat("\nemission:\n")
    print(object$model$parms.emission)
    return(invisible(object))
}


simulate.hmmspec <- function(object, nsim, seed=NULL, rand.emission=NULL,...) {

  if(!is.null(seed)) set.seed(seed)
  if(is.null(rand.emission)&is.null(object$rand.emission)) stop("rand.emission not specified")
  if(is.null(rand.emission)) rand.emission=object$rand.emission
  
  if(length(nsim)==1) {
   s1 = sim.mc(object$init,object$transition, nsim)
   x = sapply(s1, rand.emission, model=object)
    if (NCOL(x) > 1) 
        ret = list(s = s1, x = t(x), N = nsim)
    else ret = list(s = s1, x = x, N = nsim)
   class(ret) <- "hsmm.data"
   ret
  }
  else  .sim.mhmm(object,nsim,rand.emission)
}

.sim.mhmm <- function(model,N,rand.emission) {
  s1 = sim.mc(model$init,model$transition,N) #the hidden states
  x = sapply(s1, rand.emission, model)       #simulate the observed state sequence
  if(NCOL(x)>1) ret = list(s=s1,x=t(x),N=N)
  else ret = list(s=s1,x=x,N=N)

  class(ret) <- "mhsmm.data"
  ret                                            
}

hmmspec <- function(init, trans, parms.emission, dens.emission, rand.emission=NULL, mstep=NULL) {
 ans <- list(J=length(init), init = init, transition = trans, parms.emission = parms.emission,dens.emission=dens.emission,
             rand.emission=rand.emission,mstep=mstep)
 class(ans) <- "hmmspec"
 return(ans)
}

print.hmm <- function(x, ...) {
 cat("hmm object contains the following slots:\n")
 print(names(x))
 return(invisible(x))
}

.estep.hmm <- function(x,object) {
    K=nrow(object$model$transition)
    model=object$model
    N=x$N
    p = sapply(1:K,fn <- function(state) object$f(x$x,state,model))
    tmp = .C("mo_estep_hmm",a=as.double(t(model$transition)),pi=as.double(t(model$init)),p=as.double(t(p)),N=as.integer(x$N),nsequences=as.integer(length(x$N)),
      K=as.integer(K),
      alpha=double((K+1)*sum(N)) ,beta=double(K*sum(N)),gam=double(K*sum(N)),ll=double(1),PACKAGE='mhsmm')      
    list(gamma=matrix(tmp$gam,ncol=K),loglik=tmp$ll)
}

hmmfit <- function(x,start.val,mstep=mstep.norm,lock.transition=FALSE,tol=1e-08,maxit=1000) 
{
  model = start.val
  K = nrow(model$trans)
  if(mode(x)=="numeric" | mode(x)=="integer") {
    warning('x is a primitive vector.  Assuming single sequence.')
    N = NN = NROW(x)
  }
  else{
    N = NN = x$N
    x = x$x
  }
  
  if(K<2) stop("K must be larger than one.")	
  if(any(dim(model$trans)!=K)) stop("dimensions of a incorrect")
  if(length(model$init)!=K) stop("dimensions of st incorrect")
  if(NROW(x)!=sum(N)) stop("dimensions incorrect")
  if(length(N)==1) NN=1
  else NN = length(N)
  
  if(is.null(mstep)) 
    if(is.null(model$mstep)) stop("mstep not specified")
    else  mstep=model$mstep      

  f=model$dens.emission

  loglik=numeric(maxit)
  loglik=NA
  loglik[1]=-Inf  
  gam = double(K*sum(N))
  for(i in 1:maxit) {  
    p = sapply(1:K,fn <- function(state) f(x,state,model))
    if(any(apply(p,1,max)==0)) stop("Some values have 0 pdf for all states!  Check your model parameters")
    
    #estep duh	
    test = .C("mo_estep_hmm",a=as.double(t(model$transition)),pi=as.double(t(model$init)),p=as.double(t(p)),
      N=as.integer(N),nsequences=as.integer(NN), K=as.integer(K),
      alpha=double((K+1)*sum(N)) ,beta=double(K*sum(N)),gam=gam,ll=double(1),PACKAGE='mhsmm')
    #mstep
    loglik[i]=test$ll                                   
  if(i>1)    if(abs(loglik[i]-loglik[i-1])<tol) break("Converged")
#    if((loglik[i]-loglik[i-1])<(-tol)) stop(paste("loglikelihood has decreased on iteration",i))
    gam = matrix(test$gam,ncol=K)
    if(any(colSums(gam)==0)) stop("Error: at least one state has an expected number of occurences equal to 0.\n This may be caused by bad starting parameters are insufficent sample size")

    if(length(formals(mstep))==2) {
      model$parms.emission = mstep(x,gam)
    }
    else if(length(formals(mstep))==4) {
      alpha = matrix(test$alpha,ncol=K+1)
      beta = matrix(test$beta,ncol=K)
      model$parms.emission = mstep(x,gam,alpha,beta)
    }
    else {
      stop("Error: M-step function is invalid.")
    }
    if(!lock.transition) {
      model$transition=matrix(test$a,nrow=K,byrow=TRUE)
      model$init=test$pi
      model$init[model$init<0]=0
      model$transition[model$transition<0]=0
    }
  }  
  ret = list(model=model,K=K,f=f,mstep=mstep,gam=gam,loglik=loglik[!is.na(loglik)],N=N,p=gam,yhat=apply(gam,1,which.max))
  class(ret) <- "hmm"
  return(ret)	
}


predict.hmm <- function(object,newdata,method="viterbi",...) {
  if(missing(newdata)) stop("no data passed to predict!")
  else x = newdata
  if(mode(x)=="numeric" | mode(x)=="integer") {
  	warning('x is a primitive vector.  Assuming single sequence.')
  	N = NROW(x)
  	if(N<1) stop("N less than one")
  	x=list(x=x,N=NROW(x))
  }
  nseq=length(x$N) 
  	N = x$N
  	NN = cumsum(c(0,x$N))
  if(method=="viterbi") {
    nseq=1
    K = object$K
    p = sapply(1:K,fn <- function(state) object$f(x$x,state,object$model))
    p[p==0]= 1e-200
    tmp = object$model$trans
    tmp[!tmp>0] =  1e-200
    logtrans = as.double(log(t(tmp)))
    tmp = object$model$init
    tmp[!tmp>0] =  1e-20
    logpi = as.double(log(t(tmp)))
    state = integer(sum(x$N))
    loglik=0
    for(i in 1:length(x$N)) {
      tmp <- .C("viterbi_hmm",a=logtrans,pi=logpi,p=as.double(log(t(p[(NN[i]+1):NN[i+1],]))),N=as.integer(x$N[i]),NN=as.integer(nseq),K=as.integer(object$K),
                q=as.integer(rep(-1,x$N[i])),loglik=as.double(c(0)),PACKAGE='mhsmm')
      loglik=loglik+tmp$loglik
      state[(NN[i]+1):NN[i+1]] = tmp$q+1
    }
    ans <- list(s=state,x=x$x,N=x$N,loglik=loglik)
  }
  else if(method=="smoothed") {
    tmp <- .estep.hmm(x,object)
    yhat <- apply(tmp$gamma,1,which.max)
    ans <- list(s=yhat,x=x$x,N=x$N,p=tmp$gamma,loglik=tmp$loglik)
  }
  else stop("Unavailable prediction method")
  class(ans) <- "hsmm.data"
  ans
}

predict.hmmspec <- function(object,newdata,method="viterbi",...) {
  object2 <- list(model=object,K=object$J,f=object$dens.emission)
  class(object2) <- "hmm"
  predict(object2,newdata,method,...)  
}
