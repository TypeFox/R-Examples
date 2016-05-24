print.hsmmspec <- function(x, ...){
  cat("Hidden semi-Markov Model specification:\n")
  cat(sprintf("J (number of states): \n%i \n", x$J))
  cat("init:\n")
  print(x$init)
  cat ("transition matrix:\n")
  print(x$transition)
  cat("emission distribution:\n")
  print(x$parms.emission)
  cat("sojourn distribution:\n")
  print(x$sojourn)
  return(invisible(x))
}

.fitnbinom <- function(eta) {  
  shiftthresh=1e-20
  maxshift =  match(TRUE,eta>shiftthresh)
  Mtmp = tail(which(eta>shiftthresh),1)

  fun1 <- function(shift) {
    m <- weighted.mean((maxshift:Mtmp)-shift,eta[maxshift:Mtmp])
    v <- as.numeric(cov.wt(data.frame((maxshift:Mtmp)-shift),wt=eta[maxshift:Mtmp])$cov)
    size <- if (v > m) m^2/(v - m) else 100
    densfun <- function(par) sum(dnbinom((maxshift:Mtmp)-shift,size=par[1],mu=par[2],log=TRUE)*eta[maxshift:Mtmp])    
    optim(c(size,m),densfun,control=list(fnscale=-1))$value
  }
  
  shift = which.max(sapply(1:maxshift,fun1))
    m <- weighted.mean((maxshift:Mtmp)-shift,eta[maxshift:Mtmp])
    v <- as.numeric(cov.wt(data.frame((maxshift:Mtmp)-shift),wt=eta[maxshift:Mtmp])$cov)
    size <- if (v > m) m^2/(v - m) else 100
    densfun <- function(par) sum(dnbinom((maxshift:Mtmp)-shift,size=par[1],mu=par[2],log=TRUE)*eta[maxshift:Mtmp])        
    tmp = optim(c(size,m),densfun,control=list(fnscale=-1))$par
    c(shift = shift,size=tmp[1],mu=tmp[2],prob=tmp[1]/(sum(tmp)))
}


.dnbinom.hsmm.sojourn <- function(x,size,prob=NULL,shift,mu=NULL,log=FALSE) {
  if(shift<0) stop(".dnbinom.hsmm.sojourn: shift must be > 0")
  if(is.null(mu)){
    if(log) dnbinom(x-shift,size,prob,log=TRUE)
    else dnbinom(x-shift,size,prob)  
  }
  else {
    if(log) dnbinom(x-shift,size=size,mu=mu,log=TRUE)
    else dnbinom(x-shift,size=size,mu=mu)    
  }
}

.rnbinom.hsmm.sojourn <- function(n,size,prob,shift) {
  if(shift<0) stop(".dnbinom.hsmm.sojourn: shift must be > 0")
  rnbinom(n,size,prob) + shift 
}


.check.hsmmspec <- function(object) {
  if(is.null(object$dens.emission)) stop("No emission density function provided!")
  if(is.null(object$init)) stop("No initial distribution specified!")
  if(is.null(object$transition)) stop("No initial distribution specified!")
  if(is.null(object$parms.emission)) stop("No emission parameters specified!")
  if(is.null(object$sojourn)) stop("No sojourn distribution specified!")
  if(!is.null(object$sojourn$d)) if(NCOL(object$sojourn$d)!=nrow(object$transition)) stop("Inconsistent sojourn d")
  if(length(object$init)!=NROW(object$transition))    stop('length(init)!=NROW(transition)')
  if(NROW(object$transition)!=NCOL(object$transition)) stop('NROW(transition)!=NCOL(transition)')
  if(!isTRUE(all.equal(sum(diag(object$transition)),0))) stop('non-zero entry on diagonal of transition matrix')  
}

hsmmspec <- function(init,transition,parms.emission,sojourn,dens.emission,rand.emission=NULL,mstep=NULL) {
  if(is.null(dens.emission)) stop("dens.emission not specified")
  if(length(init)!=NROW(transition))    stop('length(init)!=NROW(transition)')
  if(NROW(transition)!=NCOL(transition)) stop('NROW(transition)!=NCOL(transition)')
  if(is.null(sojourn$type)) stop("Sojourn distribution type not specified.")
  if(all(sojourn$type!=c("nonparametric","gamma","poisson"))) stop(paste("Invalid sojourn type specified (",sojourn$type,")"))
  ans = list(J=length(init),init=init,transition=transition,parms.emission=parms.emission,sojourn=sojourn,rand.emission=rand.emission,dens.emission=dens.emission,mstep=mstep)
  class(ans) <- 'hsmmspec'
  .check.hsmmspec(ans)
  ans  
}

simulate.hsmmspec <- function(object, nsim, seed=NULL,rand.emission=NULL,...)
{
  right.truncate=left.truncate=0
  if(!is.null(seed)) set.seed(seed)
  if(is.null(rand.emission)&is.null(object$rand.emission)) stop("rand.emission not specified")
  if(!is.null(rand.emission)) object$rand.emission=rand.emission
  
  if(length(nsim)==1) {
    s0 = sim.mc(object$init,object$transition, nsim)
    if (object$sojourn$type == "poisson") {
#        object$d = matrix(nrow = M, ncol = object$J)
#        for (i in 1:object$J) object$d[, i] = .dpois.hsmm.sojourn(1:M, object$sojourn$lambda[i],object$sojourn$shift[i])
        fn <- function(ii) .rpois.hsmm.sojourn(1,object$sojourn$lambda[ii],object$sojourn$shift[ii])
    }
    else if (object$sojourn$type == "gamma") {
#        object$d = matrix(nrow = M, ncol = object$J)
#        for (i in 1:object$J) object$d[, i] = dgamma(1:M, object$sojourn$shape[i], 
#            scale = object$sojourn$scale[i])
        fn <- function(ii) rgamma(1,shape=object$sojourn$shape[ii],scale=object$sojourn$scale[ii])
    }
#    else if (object$sojourn$type == "lnorm" | object$sojourn$type == "lognormal") {
#        object$d = matrix(nrow = M, ncol = object$J)
#        for (i in 1:object$J) {
#          object$d[, i] = dlnorm(1:M, object$sojourn$meanlog[i],object$sojourn$s.dlog[i])
#          object$d[, i] = object$d[, i]/sum(object$d[, i])          
#        }        
#    }    
    else if (object$sojourn$type == "logarithmic") {
#        object$d = matrix(nrow = M, ncol = object$J)
#        for (i in 1:object$J) object$d[, i] = .dlog(1:M,object$sojourn$shape[i])
        fn <- function(ii) .rlog(1,object$sojourn$p[ii])        
    }
#    else if (object$sojourn$type == "nbinom") {
#        object$d = matrix(nrow = M, ncol = object$J)
#        for (i in 1:object$J) object$d[, i] = .dnbinom.hsmm.sojourn(1:M,object$sojourn$size[i],object$sojourn$prob[i],object$sojourn$shift[i])
#    }    
    else if (object$sojourn$type == "nonparametric") {
      fn <- function(ii) sample(1:nrow(object$sojourn$d), 1, prob = object$sojourn$d[,ii])
    }
    else stop("Unknown sojourn distribution")
    
    u = as.integer(round(sapply(s0, fn)))
    s1 = rep(s0, u)[(left.truncate + 1):(sum(u) - right.truncate)]
    x = sapply(s1,function(i) rand.emission(i,object))
    if (NCOL(x) > 1)
        ret = list(s = s1, x = t(x), N = NCOL(x))
    else ret = list(s = s1, x = x, N = NROW(x))
    class(ret) <- "hsmm.data"
    ret
  }
  else {           #TODO rewrite this clause
    .sim.mhsmm(nsim,object,object$sojourn$type,object$rand.emission)
  }
}

#simulate a multisequence hsmm
.sim.mhsmm <- function(N,model,sojourn.distribution="nonparametric",emission=rnorm.hsmm,left.truncate=0,right.truncate=0,M=10000)
{
  M=10000
  if(!all.equal(diag(model$transition),rep(0,model$J))) stop("Diagonal of a must be all 0")
  s0 = sim.mc(model$init,model$transition,N)#the embedded markov model
  object = model
  if(sojourn.distribution=="poisson") {
  	model$d = matrix(nrow=M,ncol=model$J)
  	for(i in 1:model$J) model$d[,i] = .dpois.hsmm.sojourn(1:M,model$sojourn$lambda[i],model$sojourn$shift[i])
  } else  if(sojourn.distribution=="gamma") {
  	model$d = matrix(nrow=M,ncol=model$J)
  	for(i in 1:model$J) model$d[,i] = dgamma(1:M,model$sojourn$shape[i],scale=model$sojourn$scale[i])
  }
  else if (sojourn.distribution=="lnorm" |sojourn.distribution=="lognormal") {
    model$d = matrix(nrow = M, ncol = model$J)
    for (i in 1:model$J) {
      model$d[, i] = dlnorm(1:M, model$sojourn$meanlog[i],model$sojourn$s.dlog[i])
      model$d[, i] = model$d[, i]/sum(model$d[, i])          
    }                
  }
  else if (sojourn.distribution == "logarithmic") {
    model$d = matrix(nrow = M, ncol = model$J)
    for (i in 1:model$J) model$d[, i] = .dlog(1:M,model$sojourn$shape[i])
  }
  else if (object$sojourn$type == "nbinom") {
    model$d = matrix(nrow = M, ncol = object$J)
    for (i in 1:object$J) model$d[, i] = .dnbinom.hsmm.sojourn(1:M,object$sojourn$size[i],object$sojourn$prob[i],object$sojourn$shift[i])
  }
  else if (object$sojourn$type == "nonparametric") {
    model$d <- model$sojourn$d
  }
  else stop("Unknown sojourn distribution")

  
  fn <- function(ii) sample(1:nrow(model$d),1,prob=model$d[,ii])
  u = sapply(s0,fn) #simulate the waiting times  
  
  s1 = rep(s0,u)[(left.truncate+1):(sum(u)-right.truncate)] #simulate actual state sequence

  NN = N
  tmp = cumsum(c(1,N))
  for(i in 1:(length(tmp)-1))
    NN[i] = sum(u[tmp[i]:(tmp[i+1]-1)])
  

  x = sapply(s1,emission,model) #simulate the observed state sequence
  if(NCOL(x)>1) ret = list(s=s1,x=t(x),N=NN)
  else ret = list(s=s1,x=x,N=NN)
  class(ret) <- "hsmm.data"
  ret
}

.rlog <- function(n,p,M=10000) sample(1:M,n,TRUE,.dlog(1:M,p))

.dlog <- function(x,p) (-1 / log(1-p)) * p^x / x

.logdistrfit <- function(wt) {
  xbar = sum(wt*(1:length(wt)))
  fn <- function(p) xbar + p/((1-p)*log(1-p))
  uniroot(fn,c(1e-10,1-1e-10))$root
}      

#estimates gamma distribution parameters using method of moments
gammafit <- function(x,wt=NULL) {
   tol = 1e-08

    if(is.null(wt)) wt = rep(1,length(x))

      tmp = cov.wt(data.frame(x),wt=wt)
      xhat = tmp$center
      xs = sqrt(tmp$cov)
      s = log(xhat) - mean(weighted.mean(log(x),wt))    
      aold = (xhat/xs)^2
      a = Inf
      while(abs(a-aold)>tol) {
        a = aold - (log(aold) - digamma(aold) - s)/((1/aold) - trigamma(aold))        
        aold=a
      }
      return(list(shape=a,scale=xhat/a))
}

sim.mc <- function(init,transition,N) {
  if(!all.equal(rowSums(transition),rep(1,nrow(transition)))) stop("Rows of a must sum to one")
  if(!all.equal(sum(init),1)) stop("st must sum to one")
	a0 =  t(apply(transition,1,cumsum))
	st= cumsum(init)
	state = integer(sum(N))
  .C("sim_mc",as.double(st),as.double(a0),as.integer(nrow(transition)),state=state,as.integer(N),as.integer(length(N)),PACKAGE='mhsmm')$state
}

#rng for a shifted poisson distrubtion
.rpois.hsmm.sojourn <- function(n,lambda,shift) rpois(n,lambda)+shift

#density function for a shifted poisson distribution
.dpois.hsmm.sojourn <- function(x=NULL,lambda,shift,log=FALSE)   {
  if(shift<0) stop(".dpois.hsmm.sojourn: shift must be > 0")
  if(log) dpois(x-shift,lambda,log=TRUE)
  else dpois(x-shift,lambda)
}

.build_d <- function(model,M) {
  sojourn.distribution=model$sojourn$type
  J = model$J
  if(sojourn.distribution=='nonparametric' | sojourn.distribution=="ksmoothed-nonparametric") {
    if(!is.null(model$sojourn$d)) {
      if(ncol(model$sojourn$d)!=J) stop("ncol(model$d)!=J")
      M = nrow(model$sojourn$d)
      model$d = model$sojourn$d
    }
    else stop("Sojourn distribution (model$sojourn$d) not specified.")
  }

  if(sojourn.distribution=="poisson") {
    if(is.null(model$sojourn$d)) {
      if(is.null(model$sojourn$lambda)) stop('Invalid waiting parameters supplied')
      if(is.null(model$sojourn$shift)) stop('No shift parameter provided for Poisson sojourn distribution (must be at least 1)')
      model$d = matrix(nrow=M,ncol=model$J)  	
      for(i in 1:J) model$d[,i] = .dpois.hsmm.sojourn(1:M,model$sojourn$lambda[i],model$sojourn$shift[i])
    }
    else
      for(i in 1:J) model$d = model$sojourn$d
  }
  
  if(sojourn.distribution=="lnorm") {
    if(is.null(model$sojourn$d)) {
      if(is.null(model$sojourn$meanlog) | is.null(model$sojourn$s.dlog)) stop('Invalid waiting parameters supplied')
      model$d = matrix(nrow=M,ncol=model$J)  	
      for(i in 1:J) model$d[,i] = dlnorm(1:M,model$sojourn$meanlog[i],model$sojourn$s.dlog[i])
    }
    else
      for(i in 1:J) model$d = model$sojourn$d
  }

  if(sojourn.distribution=="gamma") {
    if(is.null(model$sojourn$shape) | is.null(model$sojourn$scale)) {
      if(is.null(model$sojourn$d))
        stop('Invalid waiting parameters supplied')
      else model$d = model$sojourn$d
    }
    else {    
      model$d = matrix(nrow=M,ncol=model$J)
      for(i in 1:J) model$d[,i] = dgamma(1:M,shape=model$sojourn$shape[i],scale=model$sojourn$scale[i])
    }
  }  

  if(sojourn.distribution=="logarithmic") {
    if(is.null(model$sojourn$shape)) {
      if(is.null(model$sojourn$d))
        stop('Invalid waiting parameters supplied')
      else model$d = model$sojourn$d
    }
    else {    
      model$d = matrix(nrow=M,ncol=model$J)
      for(i in 1:J) model$d[,i] = .dlog(1:M,model$sojourn$shape[i])
    }
  }  
  
  if (sojourn.distribution == "nbinom") {
    model$d = matrix(nrow = M, ncol = J)
    if(is.null(model$sojourn$mu))    for (i in 1:J) model$d[, i] = .dnbinom.hsmm.sojourn(1:M,size=model$sojourn$size[i],prob=model$sojourn$prob[i],shift=model$sojourn$shift[i])
    if(is.null(model$sojourn$prob))    for (i in 1:J) model$d[, i] = .dnbinom.hsmm.sojourn(1:M,size=model$sojourn$size[i],mu=model$sojourn$mu[i],shift=model$sojourn$shift[i])
  }

  model$d=head(model$d,M)
  model$D = apply(model$d,2,function(x) rev(cumsum(rev(x))))
  model
}

hsmmfit <- function(x,model,mstep=NULL,M=NA,maxit=100,lock.transition=FALSE,lock.d=FALSE,graphical=FALSE) {
  sojourn.distribution=model$sojourn$type
  tol=1e-4
  ksmooth.thresh = 1e-20 #this is a threshold for which d(u) values to use - if we throw too many weights in the default density() seems to work quite poorly
  shiftthresh = 1e-20 #threshold for effective "0" when considering d(u)
  J = nrow(model$transition)
  model$J = J
  
  if(is.null(mstep)) 
    if(is.null(model$mstep)) stop("mstep not specified")
    else  mstep=model$mstep      

  .check.hsmmspec(model)
  
  f=model$dens.emission

  if(mode(x)=="numeric" | mode(x)=="integer") {
    warning('x is a primitive vector.  Assuming single sequence.')
    NN = NROW(x)    
  }
  else{
    NN = x$N
    x = x$x
  }
  if(is.na(M)) M = max(NN)      
  if(length(model$init)!=J) stop("length(model$init)!=J")
  if(NROW(x)!=sum(NN)) stop("NROW(x)!=sum(NN)")
  model <- .build_d(model,M)
  
  new.model = model
  ll = rep(NA,maxit)
  rm(model)

#  N.debug = list()
  for(it in 1:maxit) {
    if(graphical)   plot.hsmm(list(model=new.model,J=J))
    p = sapply(1:J,function(state) f(x,state,new.model))
    if(any(is.na(p))) stop("NAs detected in b(x), check your supplied density function")
    if(any(apply(p,1,max)==0)) stop("Some values have 0 pdf for all states!  Check your model parameters")
    
#    print(paste("Iteration",it))
# E-STEP
    B  = .C("backward",
          transition=as.double(new.model$transition),
          init=as.double(new.model$init),
          p=as.double(p),
          d=as.double(new.model$d),
          D=as.double(new.model$D),
          timelength=as.integer(NN),
          J=as.integer(J),
          M=as.integer(rep(M,J)),
          L1 = double(NROW(x)*J),N = double(NROW(x)),
          eta = double(M*J),
          F1=double(J*NROW(x)),
          si=double(J*NROW(x)),
          gamma=double(J*NROW(x)),
          nsequences=as.integer(length(NN)),
          totallength=NROW(x),
          G=double(J*NROW(x)),
          PACKAGE='mhsmm')

#M-Step

#    N.debug[[it]] = B$N
    if(any(is.nan(B$gamma))) {
      warning("NaNs detected in gamma.  Exiting...")
      return(B)
    }
    if(any(B$gamma<0)) B$gamma = zapsmall(B$gamma)      
    if(any(B$eta<0)) B$eta = zapsmall(B$eta)      
    if(any(B$N<0))  B$N = zapsmall(B$N)

#     if(any(apply(matrix(B$gamma,ncol=J),2,function(gamma) all(gamma==0)))) stop("all negative gamma detected.  Exiting...")
     old.model = new.model
#     new.model = list(parms.emission=mstep(x,matrix(B$gamma,ncol=J)))
    state_wt <- matrix(B$gamma,ncol=J)
    if(any(colSums(state_wt)==0)) stop("Error: at least one state has an expected number of occurences equal to 0.\n This may be caused by bad starting parameters are insufficent sample size")
    new.model$parms.emission = mstep(x,state_wt)

    if(lock.d) {
      new.model$d = old.model$d
      new.model$D = old.model$D
    }
    else {
      if(sojourn.distribution=="nonparametric") {
#       print("Fitting non-parametric waiting time distribution.")
       new.model$d = apply(matrix(B$eta,ncol=J),2,function(x) x/sum(x))
  
       }
       else if(sojourn.distribution=="ksmoothed-nonparametric") {
          new.model$d = apply(matrix(B$eta+1e-100,ncol=J),2,function(x) x/sum(x))    
          for(i in 1:J) {
            new.model$d[,i] = approx(density(which(new.model$d[,i]>ksmooth.thresh),weights=new.model$d[which(new.model$d[,i]>ksmooth.thresh),i],from=1,n=M),xout=1:M)$y
            new.model$d[is.na(new.model$d[,i]),i] = 0
            new.model$d[,i] = (new.model$d[,i]+1e-300)/sum(new.model$d[,i])
          }
       }
       
       else if(sojourn.distribution=="poisson") {
        #    print("Fitting Poisson waiting time distribution.")
            new.model$d = apply(matrix(B$eta,ncol=J),2,function(x) x/sum(x))
            new.model$sojourn$lambda = numeric(J)
            new.model$sojourn$shift = numeric(J)            
         for(i in 1:J) {
           eta = new.model$d[,i]
           maxshift =  match(TRUE,eta>shiftthresh)
           Mtmp = tail(which(eta>shiftthresh),1)
  	       new.model$sojourn$shift[i] = which.max(sapply(1:maxshift, function(shift) .dpois.hsmm.sojourn(x = maxshift:Mtmp,lambda=((maxshift:Mtmp)-shift)%*%eta[maxshift:Mtmp],shift=shift,log=TRUE)%*%eta[maxshift:Mtmp]))
    	     new.model$sojourn$lambda[i] = ((new.model$sojourn$shift[i]:Mtmp)-new.model$sojourn$shift[i])%*%eta[new.model$sojourn$shift[i]:Mtmp]         
    	     new.model$d[,i] = .dpois.hsmm.sojourn(1:M,new.model$sojourn$lambda[i],new.model$sojourn$shift[i])
         }
       }
       else if(sojourn.distribution=="nbinom") {       
        new.model$d = matrix(nrow=M,ncol=J)
        new.model$sojourn$size = numeric(J)
        new.model$sojourn$shift = integer(J)
        new.model$sojourn$mu = numeric(J)                            
        new.model$sojourn$prob = numeric(J)                            
        eta = matrix(B$eta,ncol=J)
        for(i in 1:J) { 
          tmp = .fitnbinom(eta[,i])          
          new.model$sojourn$shift[i] = tmp[1]
          new.model$sojourn$size[i] =  tmp[2]
          new.model$sojourn$mu[i] =  tmp[3]
          new.model$sojourn$prob[i] =  tmp[4]
    	    new.model$d[,i] =  .dnbinom.hsmm.sojourn(1:M,new.model$sojourn$size[i],new.model$sojourn$prob[i],new.model$sojourn$shift[i])
        }
       }
       else if(sojourn.distribution=="gamma") {
#            print("Fitting Gamma waiting time distribution.")
#            new.model$d = apply(matrix(B$eta+1e-100,ncol=J),2,function(x) x/sum(x))
            new.model$d = matrix(B$eta,ncol=J)
  	        new.model$sojourn$shape = numeric(J)
  	       new.model$sojourn$scale = numeric(J)
            for(i in 1:J) {           
              tmp = gammafit(1:M,wt=new.model$d[,i])
	            new.model$sojourn$shape[i] = tmp$shape
      	      new.model$sojourn$scale[i] = tmp$scale
              new.model$d[,i] = dgamma(1:M,shape=tmp$shape,scale=tmp$scale)              
            }          
       }
       else if(sojourn.distribution=="logarithmic") {
            new.model$d = apply(matrix(B$eta+1e-100,ncol=J),2,function(x) x/sum(x))
  	        new.model$sojourn$shape = numeric(J)
            for(i in 1:J) {           
    	        new.model$sojourn$shape[i] = .logdistrfit(wt=new.model$d[,i])
              new.model$d[,i] = .dlog(1:M,new.model$sojourn$shape[i])
            }
       }
       else if(sojourn.distribution=="lnorm") {
            eta = matrix(B$eta,ncol=J)
            new.model$d = matrix(nrow=M,ncol=J)
  	        new.model$sojourn$meanlog = numeric(J)
  	        new.model$sojourn$s.dlog = numeric(J)
            for(i in 1:J) {           
    	        new.model$sojourn$meanlog[i] = weighted.mean(log(1:M),eta[,i])
    	        new.model$sojourn$s.dlog[i] = sqrt(cov.wt(data.frame(log(1:M)),eta[,i])$cov)
              new.model$d[,i] = dlnorm(1:M,new.model$sojourn$meanlog[i],new.model$sojourn$s.dlog[i])
              new.model$d[,i] = new.model$d[,i]/sum(new.model$d[,i])              
            }
      }
    else stop("Invalid sojourn distribution")         
       new.model$D = apply(new.model$d,2,function(x) rev(cumsum(rev(x))))
#       new.model$D[new.model$D<1e-300]=1e-300
    }     
     if(lock.transition) {
       new.model$init=old.model$init
       new.model$transition = old.model$transition
     }
     else {
       new.model$init=B$init
       new.model$init[new.model$init<0]=0
       new.model$transition = matrix(B$transition,ncol=J)
       new.model$transition[new.model$transition<0]=0
     }
#     if(any(is.na(
     ll[it]=sum(log(B$N))
     #mstep(x,matrix(B$gamma,ncol=J),J,NROW(x))
     new.model$J = J
     if(it>2) if(abs(ll[it]-ll[it-1])<tol) break()
  }
  ## new.model$dens.emssion <- f
  ## new.model$sojourn$type = sojourn.distribution
  class(new.model) <- "hsmmspec"
  ret = list(loglik=ll[!is.na(ll)],model=new.model,B=B,M=M,J=J,NN=NN,f=f,mstep=mstep,yhat=apply(matrix(B$gamma,ncol=J),1,which.max))#,N.debug=N.debug)
  class(ret) <- "hsmm"
  ret
}

predict.hsmm <- function(object,newdata,method="viterbi",...) {
  if(missing(newdata)) stop("newdata missing!")
  else x=newdata
  J = object$J
  m=-1e300
  if(mode(x)=="numeric" | mode(x)=="integer") {
  	warning('x is a primitive vector.  Assuming single sequence.')
  	x0 = x
  	N = NROW(x)
  	NN = c(0,N)
  	if(N<1) stop("N less than one")
  }
  else{
  	N = x$N
  	NN = cumsum(c(0,x$N))
  	x0 = x$x
  }
  statehat = integer(NROW(x))
  statehat=NA
  if(method=="viterbi") {
    M = nrow(object$model$d)    
    loga = as.double(log(object$model$transition))
    loga[loga==-Inf]=m
    logstart = as.double(log(object$model$init))
    logstart[logstart==-Inf|is.nan(logstart)] = m
    d = apply(object$model$d,2,function(x) x/sum(x))
    D = apply(d,2,function(x) rev(cumsum(rev(x))))
    d = log(d)
    d[d==-Inf]=m
    D = log(D)
    D[D==-Inf]=m
    loglik=0
  for(i in 1:length(N)) {
    if(NCOL(x0)==1)    b = log(unlist(sapply(1:J,function(state) object$f(x0[(NN[i]+1):NN[i+1]],state,object$model))))
    else    b = log(unlist(sapply(1:J,function(state) object$f(x0[(NN[i]+1):NN[i+1],],state,object$model))))
    b[b==-Inf]=m
    tmp = .C("viterbi",
                a=loga,
                pi=logstart,
                p=as.double(b),
                d=as.double(d),
                D=as.double(D),
                timelength=as.integer(N[i]),
                J=as.integer(J), 
                M=as.integer(rep(M,J)),
                alpha = double(N[i]*J),
                statehat=integer(N[i]),
                psi_state0=integer(N[i]*J),
                psi_time0=integer(N[i]*J)          
                ,PACKAGE='mhsmm')
    loglik=loglik+max(tmp$alpha[N[i]*(1:3)])
    statehat[(NN[i]+1):NN[i+1]] = tmp$statehat+1
  }
  ans <- list(x=x,s=statehat,N=N,loglik=loglik)
 }
 else if(method=="smoothed") {
   M = nrow(object$model$d)    
   m <- object$model
   m$dens.emission <- object$f
   tmp <- hsmmfit(x,m,object$mstep,maxit=1,M=M)
   ans <- list(x=x$x,s=tmp$yhat,N=x$N,p=matrix(tmp$B$gamma,ncol=object$J))
  }
  else stop(paste("Unavailable prediction method",method))
    
  class(ans) <- 'hsmm.data'  
  ans
# tmp
}


.add.states <- function(states,ht=0,greyscale=FALSE,leg=NA,J=length(unique(states)),time.scale=24,shift=0) {  
  J = length(unique(states))  

  if(greyscale) cols=c("#FFFFFF" ,"#F0F0F0" ,"#D9D9D9", "#BDBDBD" ,"#969696", "#737373", "#525252", "#252525")
  else cols = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3") #kind regards to RBrewerPal for these values  

  hats = rle(states)
  hats = list(intervals=cumsum(c(0,hats$lengths))/time.scale+shift,state=hats$values)
  for(ii in 1:length(hats$state)) 
     if(greyscale)  rect(hats$intervals[ii],ht,hats$intervals[ii+1],ht+(axTicks(2)[2]-axTicks(2)[1])/5,col=cols[hats$state[ii]],border=1)
     else rect(hats$intervals[ii],ht,hats$intervals[ii+1],ht+(axTicks(2)[2]-axTicks(2)[1])/5,col=cols[hats$state[ii]],border=cols[hats$state[ii]])
  if(any(!is.na(leg))) legend("topleft",legend=leg,fill=cols,bg="white")
}


summary.hsmm <- function(object,...) {
	cat("\nStarting distribution = \n")
	print(object$model$init,2)	
	cat("\nTransition matrix = \n")
	print(object$model$transition,2)
	cat("\nSojourn distribution parameters = \n")
	print(object$model$sojourn)
	cat("\nEmission distribution parameters = \n")
	print(object$model$parms.emission)
}


predict.hsmmspec <- function(object,newdata,method="viterbi",M=NA,...) {
  if(class(newdata)=="hsmm.data") NN = newdata$N
  else NN = length(newdata)
  if(is.na(M)) M = max(NN)
  .check.hsmmspec(object)
  model <- .build_d(object,M)
  object2 <- list(model=model,J=model$J,f=model$dens.emission)
  class(object2) <- "hsmm"
  predict(object2,newdata,method)    
}

