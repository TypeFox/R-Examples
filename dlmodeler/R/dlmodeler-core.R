# Author: cns
# core functions

dlmodeler.build <- 
function(a0=NULL, P0=NULL, P0inf=NULL,
         Tt=NULL, Rt=NULL, Qt=NULL,
         Zt=NULL, Ht=NULL,
         dimensions=NULL,
         name='noname', components=list())
{
	if( is.null(dimensions) )
	{
		components[[name]] <- matrix(1,NROW(Tt),1)
		mod <- list(
			a0=a0, P0=P0, P0inf=P0inf,
			Tt=Tt, Rt=Rt, Qt=Qt,
			Zt=Zt, Ht=Ht,
			name=name,
			components=components)
	}
	else
	{
		if( length(dimensions)!=3 ) stop("dimensions are invalid, they should have the form: c(m,r,d)")
		m <- dimensions[1] # dim state vector
		r <- dimensions[2] # dim state disturbances
		d <- dimensions[3] # dim observation vector
		components[[name]] <- matrix(1,m,1)
		mod <- list(
			a0=matrix(0,m,1),
			P0=matrix(0,m,m),
			P0inf=matrix(0,m,m),
			Tt=matrix(0,m,m),
			Rt=matrix(0,m,r),
			Qt=matrix(0,r,r),
			Zt=matrix(0,d,m),
			Ht=matrix(0,d,d),
			name=name,
			components=components)
	}
	class(mod) <- 'dlmodeler' # a nice class name!
	return(mod)
}



dlmodeler.check <-
function(model, yt=NULL)
{
	if( class(model)!='dlmodeler' ) stop("model should be of class 'dlmodeler'")
	
	t <- FALSE        # time varying model
	n <- rep(NA,5)    # number of time steps (in case of time-varying model)
	m <- NROW(model$Tt) # dim of state vector
	r <- NROW(model$Qt) # dim of state disturbance
	d <- NROW(model$Zt) # dim of observation vector
	s <- TRUE         # status (TRUE if model is valid)
	
	# check the dimensions of the initial conditions
	if( NROW(model$a0)!=m ) { warning("Inconsistent model: dimension of a0 should be ",m," x 1"); s <- FALSE
	} else if( NCOL(model$a0)!=1 ) { warning("Inconsistent model: dimension of a0 should be ",m," x 1"); s <- FALSE }
	if( NROW(model$P0)!=m ) { warning("Inconsistent model: dimension of P0 should be: ",m," x ",m); s <- FALSE
	} else if( NCOL(model$P0)!=m ) { warning("Inconsistent model: dimension of P0 should be: ",m," x ",m); s <- FALSE }
	if( NROW(model$P0inf)!=m ) { warning("Inconsistent model: dimension of P0inf should be: ",m," x ",m); s <- FALSE
	} else if( NCOL(model$P0inf)!=m ) { warning("Inconsistent model: dimension of P0inf should be: ",m," x ",m); s <- FALSE }
	
	# check the dimensions of the state equation
	if( NROW(model$Tt)!=m ) { warning("Inconsistent model: dimension of Tt should be: ",m," x ",m); s <- FALSE 
	} else if( NCOL(model$Tt)!=m ) { warning("Inconsistent model: dimension of Tt should be: ",m," x ",m); s <- FALSE }
	if( NROW(model$Rt)!=m ) { warning("Inconsistent model: dimension of Rt should be: ",m," x ",r); s <- FALSE 
	} else if( NCOL(model$Rt)!=r ) { warning("Inconsistent model: dimension of Rt should be: ",m," x ",r); s <- FALSE }
	if( NROW(model$Qt)!=r ) { warning("Inconsistent model: dimension of Qt should be: ",r," x ",r); s <- FALSE 
	} else if( NCOL(model$Qt)!=r ) { warning("Inconsistent model: dimension of Qt should be: ",r," x ",r); s <- FALSE }	
	
	# check the dimensions of the observation equation
	if( NROW(model$Zt)!=d ) { warning("Inconsistent model: dimension of Zt should be: ",d," x ",m); s <- FALSE 
	} else if( NCOL(model$Zt)!=m ) { warning("Inconsistent model: dimension of Zt should be: ",d," x ",m); s <- FALSE }
	if( NROW(model$Ht)!=d ) { warning("Inconsistent model: dimension of Ht should be: ",d," x ",d); s <- FALSE 
	} else if( NCOL(model$Ht)!=d ) { warning("Inconsistent model: dimension of Ht should be: ",d," x ",d); s <- FALSE }
	
	# check the time-varying components
	if( length(dim(model$Tt))==3 ) { n[1] <- dim(model$Tt)[3]; t <- TRUE }
	if( length(dim(model$Rt))==3 ) { n[2] <- dim(model$Rt)[3]; t <- TRUE }
	if( length(dim(model$Qt))==3 ) { n[3] <- dim(model$Qt)[3]; t <- TRUE }
	if( length(dim(model$Zt))==3 ) { n[4] <- dim(model$Zt)[3]; t <- TRUE }
	if( length(dim(model$Ht))==3 ) { n[5] <- dim(model$Ht)[3]; t <- TRUE }
	if( t ) if (min(n,na.rm=TRUE)!=max(n,na.rm=TRUE)) { warning("Inconsistent model: dimension of time-varying terms (Tt, Rt, Qt, Zt, Ht) should match: ",n); s <- FALSE }
	
	# check the components
	for( compname in names(model$components) ) {
		comp <- model$components[[compname]]
		if( NROW(comp)!=m ) { warning("Inconsistent model: dimension of component ",compname," should be ",m," x 1"); s <- FALSE
		} else if( NCOL(comp)!=1 ) { warning("Inconsistent model: dimension of component ",compname," should be ",m," x 1"); s <- FALSE }
		if( apply(comp==1|comp==0,2,"sum")!=m ) { warning("Inconsistent model: ",compname," be a vector of zeros and ones"); s <- FALSE }
	}
	if( apply(model$components[[model$name]]==1,2,"sum")!=m ) { warning("Inconsistent model: component named ",model$name," should be a vector of ones"); s <- FALSE }
	
	# check compatibility with data
	if( !is.null(yt) ) {
		if( !is.matrix(yt) ) yt <- matrix(yt,nrow=1)
		if( NROW(yt)!=d ) { warning("Inconsistent data: dimension of yt should be: ",d," x ",'n'); s <- FALSE
		} else if( t & NCOL(yt)!=max(n,na.rm=TRUE) ) { warning("Inconsistent data: dimension of yt should be: ",d,"x",max(n,na.rm=TRUE)); s <- FALSE }
	}
	
	return(list(status=s, m=m, r=r, d=d, timevar=t,
				timevar.Tt=n[1], timevar.Rt=n[2], timevar.Qt=n[3],
				timevar.Zt=n[4], timevar.Ht=n[5]))
}



print.dlmodeler <-
function(x,...)
{
	if( class(x)!='dlmodeler' ) stop("argument should be of class 'dlmodeler'")
	
	mdl.dim <- dlmodeler.check(x)
	if( mdl.dim$timevar==FALSE ) {
		cat("constant dlmodel(state dim=",mdl.dim$m,", dist dim=",mdl.dim$r,", obs dim=",mdl.dim$d,") '",x$name,"'\n",sep='')
	} else {
		n <- max(mdl.dim$timevar.Tt,mdl.dim$timevar.Rt,mdl.dim$timevar.Qt,mdl.dim$timevar.Zt,mdl.dim$timevar.Ht,na.rm=TRUE)
		cat("time-varying dlmodel(state dim=",mdl.dim$m,", dist dim=",mdl.dim$r,", obs dim=",mdl.dim$d,", time steps=",n,") '",x$name,"'\n",sep='')
		if( !is.na(mdl.dim$timevar.Tt) ) cat(" - state transition matrix is time-variant\n")
		if( !is.na(mdl.dim$timevar.Rt) ) cat(" - state disturbance selection matrix is time-variant\n")
		if( !is.na(mdl.dim$timevar.Qt) ) cat(" - state disturbance covariance matrix is time-variant\n")
		if( !is.na(mdl.dim$timevar.Zt) ) cat(" - observation design matrix is time-variant\n")
		if( !is.na(mdl.dim$timevar.Ht) ) cat(" - observation disturbance covariance matrix is time-variant\n")
	}
	
	nb.components <- length(x$components)
	if( nb.components>1 ) {
		cat(" - model has",nb.components,"components: ")
		cat(names(x$components),sep=", ")
		cat("\n")
	} else if( nb.components==1 ) {
		cat(" - model has",nb.components,"component: ")
		cat(names(x$components),sep=", ")
		cat("\n")
	}
	
	x <- dlmodeler.build.unknowns(x)
	if( length(x$unknowns)==1 ) {
		cat(" - model has",length(x$unknowns),"unknown parameter\n")
	} else if( length(x$unknowns)>1 ) {
		cat(" - model has",length(x$unknowns),"unknown parameters\n")
	}
	
	invisible(x)
}



dlmodeler.timevar.fun <-
function(x, y, fun)
{
	ndim.x <- length(dim(x))
	ndim.y <- length(dim(y))
	if( ndim.x==1 & ndim.y==1 ) return(fun(x,y)) # ok easy one
	if( ndim.x==2 & ndim.y==2 ) return(fun(x,y)) # also trivial
	if( ndim.x==3 & ndim.y==2 ) {
		# a bit trickier: only one parameter is time-varying
		# deducing the resulting size from the first element
		N <- dim(x)[3]
		# for some strange reason, need to coerce the sub matrix
		# into a matrix...
		guess <- fun(matrix(x[,,1],NROW(x),NCOL(x)),y)
		ret <- array(dim=c(NROW(guess),NCOL(guess),N))
		for( i in 1:N ) ret[,,i] <- fun(matrix(x[,,i],NROW(x),NCOL(x)),y)
		return(ret)
	}
	if( ndim.x==2 & ndim.y==3 ) {
		# same as above
		N <- dim(y)[3]
		guess <- fun(x,matrix(y[,,1],NROW(y),NCOL(y)))
		ret <- array(dim=c(NROW(guess),NCOL(guess),N))
		for( i in 1:N ) ret[,,i] <- fun(x,matrix(y[,,i],NROW(y),NCOL(y)))
		return(ret)
	}
	if( ndim.x==3 & ndim.y==3 ) {
		# same as above, except the two parameters are now time varying as well
		N <- dim(x)[3] # assuming = dim(y)[3]
		guess <- fun(matrix(x[,,1],NROW(x),NCOL(x)),matrix(y[,,1],NROW(y),NCOL(y)))
		ret <- array(dim=c(NROW(guess),NCOL(guess),N))
		for( i in 1:N ) ret[,,i] <- fun(matrix(x[,,i],NROW(x),NCOL(x)),matrix(y[,,i],NROW(y),NCOL(y)))
		return(ret)
	}
	stop("Don't know how to handle dimensions: ",ndim.x,"x",ndim.y)
}



dlmodeler.bdiag <-
function(x, y)
{
	# here is graphically what this function does:
	rbind( cbind(             x             , matrix(0,NROW(x),NCOL(y)) ),
	       cbind( matrix(0,NROW(y),NCOL(x)) ,            y              )
	)
}



dlmodeler.add <-
function(e1, e2, name=NULL)
{
	if( is.null(e1) ) return(e2)
	if( is.null(e2) ) return(e1)
  if( class(e1)=='numeric' | class(e1)=='integer' )
    e1 <- dlmodeler.build.constant(e1)
	if( class(e2)=='numeric' | class(e2)=='integer' )
	  e2 <- dlmodeler.build.constant(e2)
	if( class(e1)!='dlmodeler' ) stop("e1 should be of class 'dlmodeler'")
	if( class(e2)!='dlmodeler' ) stop("e2 should be of class 'dlmodeler'")
	if( is.null(name) ) name <- paste(e1$name,e2$name,sep='+')
	
	# concatenate the state vectors
	a0 <- rbind(e1$a0,e2$a0)
	P0 <- dlmodeler.bdiag(e1$P0,e2$P0)
	P0inf <- dlmodeler.bdiag(e1$P0inf,e2$P0inf)
	Tt <- dlmodeler.timevar.fun(e1$Tt,e2$Tt,dlmodeler.bdiag)
	Rt <- dlmodeler.timevar.fun(e1$Rt,e2$Rt,dlmodeler.bdiag)
	Qt <- dlmodeler.timevar.fun(e1$Qt,e2$Qt,dlmodeler.bdiag)
	
	# add outputs together
	Zt <- dlmodeler.timevar.fun(e1$Zt,e2$Zt,cbind)
	myaddition <- function(x,y) (x+y)
	Ht <- dlmodeler.timevar.fun(e1$Ht,e2$Ht,myaddition)
	
	# keep track of components
	comp1 <- lapply(e1$components,function(x) rbind(x,matrix(0,NROW(e2$Tt),1)) )
	comp2 <- lapply(e2$components,function(x) rbind(matrix(0,NROW(e1$Tt),1),x) )
	components <- c(comp1,comp2)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name,components=components))
}

"+.dlmodeler" <- function(e1, e2) dlmodeler.add(e1, e2)



dlmodeler.bind <-
function(e1, e2, name=NULL)
{
	if( is.null(e1) ) return(e2)
	if( is.null(e2) ) return(e1)
	if( class(e1)!='dlmodeler' ) stop("e1 should be of class 'dlmodeler'")
	if( class(e2)!='dlmodeler' ) stop("e2 should be of class 'dlmodeler'")
	if( is.null(name) ) name <- paste(e1$name,e2$name,sep='&')
	
	# concatenate the state vectors
	a0 <- rbind(e1$a0,e2$a0)
	P0 <- dlmodeler.bdiag(e1$P0,e2$P0)
	P0inf <- dlmodeler.bdiag(e1$P0inf,e2$P0inf)
	Tt <- dlmodeler.timevar.fun(e1$Tt,e2$Tt,dlmodeler.bdiag)
	Rt <- dlmodeler.timevar.fun(e1$Rt,e2$Rt,dlmodeler.bdiag)
	Qt <- dlmodeler.timevar.fun(e1$Qt,e2$Qt,dlmodeler.bdiag)
	
	# concatenate the outputs together
	Zt <- dlmodeler.timevar.fun(e1$Zt,e2$Zt,dlmodeler.bdiag)
	Ht <- dlmodeler.timevar.fun(e1$Ht,e2$Ht,dlmodeler.bdiag)
	
	# keep track of components
	comp1 <- lapply(e1$components,function(x) rbind(x,matrix(0,NROW(e2$Tt),1)) )
	comp2 <- lapply(e2$components,function(x) rbind(matrix(0,NROW(e1$Tt),1),x) )
	components <- c(comp1,comp2)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name,components=components))
}

"%%.dlmodeler" <- function(e1, e2) dlmodeler.bind(e1, e2)



dlmodeler.multiply <-
function(e1, e2)
{
  if( is.null(e1) ) return(e2)
  if( is.null(e2) ) return(e1)
  
  if( class(e1)=='dlmodeler' ) {
    mod <- e1
    sca <- e2
  } else if( class(e2)=='dlmodeler' ) {
    mod <- e2
    sca <- e1
  } else stop("either e1 or e2 should be of class 'dlmodeler'")
  
  if( class(sca)!='numeric' & class(sca)!='integer' )
    stop("operand should be of class 'numeric'")
  
  mod$Zt <- mod$Zt * sca
  
  return(mod)
}    

"*.dlmodeler" <- function(e1, e2) dlmodeler.multiply(e1, e2)



dlmodeler.filter <-
function(yt, model, backend=c('KFAS','FKF','dlm'), smooth=FALSE, raw.result=FALSE, logLik=FALSE, filter=TRUE)
{
	if( !is.matrix(yt) ) yt <- matrix(yt,nrow=1)
	if( class(model)!='dlmodeler' ) stop("model should be of class 'dlmodeler'")
	if( backend[1]=='KFAS') { f <- dlmodeler.filter.KFAS(yt,model,raw.result|smooth,logLik,filter)
	} else if( backend[1]=='FKF') { f <- dlmodeler.filter.FKF(yt,model,raw.result|smooth,logLik,filter)
	} else if( backend[1]=='dlm') { f <- dlmodeler.filter.dlm(yt,model,raw.result|smooth,logLik,filter)
	} else stop("Unknown backend: ",backend[1])
	class(f) <- 'dlmodeler.filtered'
	if( smooth ) {
		f$smoothed <- dlmodeler.smooth(f,raw.result)
	}
	return(f)
}



dlmodeler.smooth <-
function(filt, raw.result=FALSE)
{
	if( class(filt)!='dlmodeler.filtered' ) stop("argument should be of class 'dlmodeler.filtered'")
	if( is.null(filt$raw.result) ) stop("Raw results are not available; use raw.result=TRUE when calling dlmodeler.filter()")
	if( filt$backend=='KFAS') { s <- dlmodeler.smooth.KFAS(filt)
	} else if( filt$backend=='FKF') { s <- dlmodeler.smooth.FKF(filt)
	} else if( filt$backend=='dlm') { s <- dlmodeler.smooth.dlm(filt)
	} else stop("Unknown backend: ",filt$backend)
	class(s) <- 'dlmodeler.smoothed'
	return(s)
}



dlmodeler.extract <-
		function(fs, model, compnames=NULL, type=c("observation","state"), value=c("mean","covariance","interval"), prob=.90)
{
	if( class(fs)!='dlmodeler.filtered' & class(fs)!='dlmodeler.smoothed' | is.null(fs$at) ) stop("Argument must be a result from dlmodeler.filter() or dlmodeler.smooth()")
	if( class(model)!='dlmodeler' | is.null(model$Zt) ) stop("model should be of class 'dlmodeler'")
	if( type[1]!='observation' & type[1]!='state' ) stop("unknown type: ",type[1])
	if( value[1]!='mean' & value[1]!='covariance' & value[1]!='interval' ) stop("unknown value: ",value[1])
	
	if( value[1]=='interval' ) {
		# extract prediction intervals...
		if( prob<=0 | prob>=1 ) stop("prob should be in the range (0 ; 1)")
		sdfact <- qnorm(.5+prob/2)
		
		# this could be optimized
		the.mean <- dlmodeler.extract(fs,model,compnames,type,value="mean")
		the.sd <- dlmodeler.extract(fs,model,compnames,type,value="covariance")
		fun.sd <- function(x) sqrt(apply(x,c(3),"diag"))
		the.sd <- lapply(the.sd,fun.sd)
		ret <- list()
		for(compname in names(the.mean)) {
			ret[[compname]] <- list(
					lower=the.mean[[compname]]-sdfact*the.sd[[compname]],
					mean=the.mean[[compname]],
					upper=the.mean[[compname]]+sdfact*the.sd[[compname]])
		}
		return(ret)
	} else {
		# which components to extract? if the list is empty, extract
		# everything and give it the model name because we always want
		# this function to return something useful
		if( is.null(compnames) ) components <- model$components else components <- model$components[compnames]
		if( length(components)==0 ) components[model$name] <- matrix(1,NROW(model$Zt),NCOL(model$Zt))
		
		# compute the means and covariances
		# state mean        : E[alpha(t)]   = a'(t)
		# output mean       : E[y(t)]       = Z'(t) %*% a'(t)
		# state covariance  : cov[alpha(t)] = P'(t)
		# output covariance : cov[y(t)]     = Z'(t) %*% P'(t) %*% t(Z'(t)) + H'(t)
		fun <- function(cmp) {
			comps <- which(cmp!=0)
			nb.comps <- length(comps)
			n <- NCOL(fs$at) # assuming = dim(fs$Pt)[3]
			
			if( value[1]=='mean' ) {
				# aprim is needed to compute state & output means
				aprim <- matrix(fs$at[comps,],nb.comps,n)
				if( type[1]=='state') return(aprim)
			} else {
				# Pprim is needed to compute state & output covariances
				Pprim <- array(fs$Pt[comps,comps,],dim=c(nb.comps,nb.comps,n))
				if( type[1]=='state') return(Pprim)
			}
			
			# the remainder of this function handles computations for outputs
			# Zprim is needed to compute both output mean & covariance
			d <- NROW(model$Zt)
			if( length(dim(model$Zt))==3 ) {
				Zprim <- array(model$Zt[,comps,],dim=c(d,nb.comps,n))
			} else {
				Zprim <- matrix(model$Zt[,comps],d,nb.comps)
			}
			if( value[1]=='mean' ) {
				# slight modification to aprim
				aprim <- array(aprim,dim=c(nb.comps,1,NCOL(fs$at)))
				ret <- dlmodeler.timevar.fun(Zprim,aprim,match.fun("%*%"))
				return( matrix(ret,d,n) )
			} else {
				fun.zpz <- function(z,p) z %*% p %*% t(z)
				ret <- dlmodeler.timevar.fun(Zprim,Pprim,fun.zpz)
				ret <- dlmodeler.timevar.fun(ret,model$Ht,match.fun("+"))
				return( ret )
			}
		}
		return(lapply(components,fun)) # lapply is magic!
	}
}




dlmodeler.forecast <-
		function(yt, model, ahead=1, iters=1, step=1, start=1, prob=.90, backend=c('KFAS','FKF','dlm'), debug=FALSE)
{
	if( !is.matrix(yt) ) yt <- matrix(yt,nrow=1)
	if(NROW(yt)>1) stop("Multivariate case not implemented yet") # TODO
	if(dlmodeler.check(model)$timevar==TRUE) stop("Time-varying models are not implemented yet") # TODO
	if( class(model)!='dlmodeler' ) stop("model should be of class 'dlmodeler'")
	if( ahead<1 ) stop("ahead should be >= 1")
	if( iters<1 ) stop("iters should be >= 1")
	if( step<1 ) stop("step should be >= 1")
	if( start<1 ) stop("start should be >= 1")
	if( prob<=0 | prob>=1 ) stop("prob should be in the range (0 ; 1)")
	
	sdfact <- qnorm(.5+prob/2)
	mymodel <- model
	pos <- start # first observation to forecast
	nas <- matrix(NA,1,ahead) # number of observations to forecast
	stp <- matrix(0,1,step) # number of observations to forecast
	simus <- data.frame()
	
	for(i in 1:iters) {
		if( i==1 | debug ) {
			# sloooow debug code
			if(pos>1) {
				obs <- matrix(yt[1,1:(pos-1)],1,pos-1)
				obs <- cbind(obs,nas)
			} else {
				obs <- cbind(nas,step)
			}
			## cat('a0:', mymodel$a0,'\n')
			## cat('P0:', mymodel$P0,'\n')
			## cat('obs:', obs,'\n')
			
			modfilter <- dlmodeler.filter(obs,mymodel,backend)
			modfilter.std <- dlmodeler.extract(modfilter,mymodel,mymodel$name,type="observation",value="covariance")
			modfilter.std <- modfilter.std[[mymodel$name]]
			modfilter.std <- matrix(sqrt(modfilter.std),nrow=1)
			
			## cat('f:', modfilter$f,'\n')
			## cat('last a:', modfilter$at,'\n')
			## cat('last P:', modfilter$Pt,'\n')
			
			sim <- data.frame(
					index=(pos):(pos+ahead-1),
					distance=1:(ahead),
					lower=t(modfilter$f)[(pos):(pos+ahead-1)]-sdfact*t(modfilter.std)[(pos):(pos+ahead-1)],
					yhat=t(modfilter$f)[(pos):(pos+ahead-1)],
					upper=t(modfilter$f)[(pos):(pos+ahead-1)]+sdfact*t(modfilter.std)[(pos):(pos+ahead-1)],
					y=yt[(pos):(pos+ahead-1)])
			simus <- rbind(simus,sim)
		} else {
			# fast restart code
			obs <- matrix(yt[1,(pos-step):(pos-1)],1,step) # assuming pos-step >= 1
			obs <- cbind(obs,nas)
			
			# update a0, P0 and P0inf
			if( i==2 ) {
				last <- pos-step
				if( pos<=modfilter$d ) stop("Exact diffuse initialization algorithm must not be interrupted, either increase the value of 'start', run this function with 'debug=TRUE' or disable exact diffuse initialization")
				mymodel$P0inf <- matrix(0,NROW(mymodel$P0inf),NCOL(mymodel$P0inf))
			} else last <- step+1
			mymodel$a0 <- modfilter$at[,last]
			mymodel$P0 <- modfilter$Pt[,,last]
			
			## cat('last i:', last,'\n')
			## cat('a0:', mymodel$a0,'\n')
			## cat('P0:', mymodel$P0,'\n')
			## cat('obs:', obs,'\n')
			
			modfilter <- dlmodeler.filter(obs,mymodel,backend)
			modfilter.std <- dlmodeler.extract(modfilter,mymodel,mymodel$name,type="observation",value="covariance")
			modfilter.std <- modfilter.std[[mymodel$name]]
			modfilter.std <- matrix(sqrt(modfilter.std),nrow=1)
			
			## cat('f:', modfilter$f,'\n')
			## cat('last a:', modfilter$at,'\n')
			## cat('last P:', modfilter$Pt,'\n')
			
			sim <- data.frame(
					index=(pos):(pos+ahead-1),
					distance=1:(ahead),
					lower=t(modfilter$f)[2:(ahead+1)]-sdfact*t(modfilter.std)[2:(ahead+1)],
					yhat=t(modfilter$f)[2:(ahead+1)],
					upper=t(modfilter$f)[2:(ahead+1)]+sdfact*t(modfilter.std)[2:(ahead+1)],
					y=yt[(pos):(pos+ahead-1)])
			simus <- rbind(simus,sim)
		}
		## print(sim)
		pos <- pos+step
	}
	return(simus)
}




