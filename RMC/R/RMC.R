# This is package RMC 

"calc.estat.vec" <-
function(probs, lambda) 
{
  nas <- is.na( probs[,1])
  
  res <- rep( NA, nrow( probs))
  res[!nas] <- apply( probs[!nas,], 1, function(x){ sum( x * ( ( 1/x)^lambda - 1))})
 
  return( res)

}


"calc.ustat.vec" <-
function(states, probs, lambda) 
{
  nas <- is.na( probs[,1])

  ref.mat <- matrix( c( ( 1:length( states))[!nas], states[!nas]), ncol=2)
  rel.probs <- probs[ref.mat]
  new.probs <- rep( NA, length( states))
  new.probs[!nas] <- rel.probs

  res <- ( 1/rel.probs)^lambda - 1

  return( res)

}


"dataEG1.doc" <-
function () 
{
}



"dataEG2.doc" <-
function () 
{
}



"dataEG3movement.doc" <-
function () 
{
}



"dataEG3patch.doc" <-
function () 
{
}



"dataEG4.doc" <-
function () 
{
}



"diagnos" <-
function( obs.states, chain.id, X=NULL, fit) 
{
	chain.states <- obs.states
	if( is.null( X))
		X <- matrix( 1, nrow=length( chain.states), ncol=1)
		
	doDist <- 1

	chain.states <- chain.states[fit$stuff$id.vec!=1]
	chain.states <- ifelse( is.na( chain.states), -9, chain.states)
	chain.id <- chain.id[fit$stuff$id.vec!=1]
	X <- data.matrix( X[fit$stuff$id.vec!=1,,drop=FALSE])

	include.ends <- FALSE
	
	temp <- .C( "diagnostics", as.integer( chain.states), as.integer( fit$stuff$n.obs), as.double( X), as.integer( fit$stuff$n.cats), as.integer( fit$stuff$n.covars), as.integer( length( chain.states)), as.integer( length( fit$stuff$n.obs)), as.double( fit$pars), as.integer( fit$stuff$phi.id), as.integer( fit$stuff$pi.id), as.integer( include.ends), jump=as.integer( rep( -99, fit$stuff$n.cats^2)), patch=as.double( rep( -99, fit$stuff$n.cats*length( chain.states))), len=integer( fit$stuff$n.cats), RjumpDists=as.double( rep( -99., sum( fit$stuff$n.obs) * fit$stuff$n.cats + 2*(fit$stuff$n.cats^2))), RdistCuts=as.integer( rep( -99, fit$stuff$n.cats^2)), as.integer( doDist),PACKAGE="RMC")
	
	res <- list()
	res[["patch"]] <- list()
	place <- 0
	for( ii in 1:fit$stuff$n.cats){
		if ( temp$len[ii]>0){
			res[["patch"]][[ii]] <- temp$patch[place + 1:temp$len[ii]]
			res[["patch"]][[ii]] <- res[["patch"]][[ii]][!is.na( res[["patch"]][[ii]])]
		}
		else
			res[["patch"]][[ii]] <- numeric(0)
		place <- place + temp$len[ii]
	}
	
	res[["jumps"]] <- matrix( temp$jump, nrow=fit$stuff$n.cats, ncol=fit$stuff$n.cats, byrow=TRUE)
	kount <- 1
	eles <- rep( -9, fit$stuff$n.cats^2)
	if( doDist==1){
		place.start <- place.end <- 0
		res[["jumpDists"]] <- list()
		res[["jumpURV"]] <- matrix( -9, nrow=fit$stuff$n.cats, ncol=fit$stuff$n.cats)
		for( ii in 1:fit$stuff$n.cats)
			res[["jumpDists"]][[ii]] <- list()
		for( ii in 1:fit$stuff$n.cats)
			for( jj in 1:fit$stuff$n.cats){
				eles[kount] <- ( jj-1) * fit$stuff$n.cats + ii
				place.end <- sum( temp$RdistCuts[1:eles[kount]])	#sum( temp$RdistCuts[eles[eles!=-9]])		#sum( temp$RdistCuts[1:ele])
				place.start <- place.end - temp$RdistCuts[eles[kount]] + 1
				kount <- kount +1
				res[["jumpDists"]][[jj]][[ii]] <- temp$RjumpDists[place.start:place.end][-1]
				if( res[["jumps"]][jj,ii] != 0)
					res[["jumpURV"]][jj,ii] <- sum( res[["jumpDists"]][[jj]][[ii]][1:res[["jumps"]][jj,ii]]) + 0.5*res[["jumpDists"]][[jj]][[ii]][res[["jumps"]][jj,ii]+1]
				else
					res[["jumpURV"]][jj,ii] <- 0.5*res[["jumpDists"]][[jj]][[ii]][1]
			}	
	}
	else
		res[["jumpDists"]] <- NULL	

	temp <- list()
	temp[["patch"]] <- res[["patch"]]
	temp[["movement"]] <- res[["jumpURV"]]
	temp[["njumps"]] <- res[["jumps"]]
	return( temp)
}



"diagnos.envel" <-
function ( obs.states, chain.id, X=NULL, fit, perc=c( 0.025, 0.975), B=100, contr=list( print.iter=50))
{
	obs.states <- as.numeric( obs.states)
	chain.id <- as.numeric( chain.id)
	
	r.chain.id <- chain.id[ fit$stuff$id.vec!=1]
	r.obs.states <- obs.states[ fit$stuff$id.vec!=1]

	n.chains <- fit$stuff$n.chains
	n.obs <- summary( as.factor( r.chain.id))
	n.cats <- fit$stuff$n.cats
	n.covars <- fit$stuff$n.covars
	
	gmat <- rep( 0, n.cats*n.covars)
	gmat[fit$stuff$phi.id==1] <- fit$par[1:sum( fit$stuff$phi.id)]
	gmat <- matrix( gmat, nrow=n.covars, ncol=n.cats)
	bmat <- matrix( rep( 0, n.cats*n.covars), nrow=n.covars, ncol=n.cats)
	bmat[fit$stuff$pi.id==1,-1] <- fit$pars[-(1:sum( fit$stuff$phi.id))]

	if( is.null( X))
		X <- matrix( 1, nrow=length( chain.id), ncol=1)
	else
		X <- data.matrix( X)
 	r.X <- data.matrix( X[ fit$stuff$id.vec!=1,])
	 	
	res <- array( -99, c( n.cats, length( r.chain.id), B), list( paste( "cat", 1:n.cats, sep=""), paste( "obs",1:length( r.chain.id),sep=""), paste( "sim", 1:B, sep="")))
	quants <- array( -99, c( n.cats, length( r.chain.id), B), list( paste( "cat", 1:n.cats, sep=""), paste( "inter",1:length( r.chain.id),sep=""), paste( "sim", 1:B, sep="")))
	cross <- array( -99, c( B, n.cats, n.cats))
	jumpDistRes <- matrix( -99, nrow=B, ncol=n.cats^2-n.cats)
	
	for( ii in 1:B){
		if( as.integer( ii/contr$print.iter) == round( ii/contr$print.iter, 8))
			print( c( "Simulation", ii))
		dat <- sim.chain( n.chains = n.chains, n.obs = n.obs, n.cats = n.cats, n.covars = n.covars, beta = bmat, gamma = gmat, X = r.X)
		dat[ is.na( r.obs.states), 2] <- NA
		temp <- diagnos( obs.states=dat[,2], chain.id=chain.id, X=X, fit=fit)
		cross[ ii,,] <- temp[["njumps"]]
		extraTemp <- c( temp$movement[ upper.tri( temp$movement, diag=F)], temp$movement[ lower.tri( temp$movement, diag=F)])
		extraTemp <- ifelse( is.na( extraTemp), 99, extraTemp)
		jumpDistRes[ii,] <- sort( extraTemp)
		for( kk in 1:n.cats)
			if( length( temp[["patch"]][[kk]])!=0){
				res[kk,1:length( temp[["patch"]][[kk]]),ii] <- sort( temp[["patch"]][[kk]])
				quants[kk,1:length( temp[["patch"]][[kk]]),ii] <- ( 1:length( temp[["patch"]][[kk]]) - 0.5) / length( temp[["patch"]][[kk]])
			}
	}
	print( "finished bootstrapping")
	jumpDistRes[jumpDistRes==99] <- NA	#careful about double operations on Matrix

	obs.states[ is.na( obs.states)] <- -9
	temp <- diagnos( obs.states, chain.id, X=X, fit)
	print( "observed residuals calculated")

	envel <- list()
	for( ii in 1:n.cats){
		tres <- as.double( res[ii,,])
		tres <- tres[tres!=-99]
		tquants <- as.double( quants[ii,,])
		tquants <- tquants[tquants!=-99]
		le <- length( temp[["patch"]][[ii]])
		if( le != 0){
			envel[[ii]] <- matrix( -99, nrow=le, ncol=4)
			colnames( envel[[ii]]) <- c( "unif", "obs", "low", "uppe")
			envel[[ii]][,"obs"] <- sort( temp[["patch"]][[ii]])
			envel[[ii]][,"unif"] <- (1:le-0.5)/le
			for( jj in 1:le){
				intev <- c( ( jj-1)/le, jj/le)
				envel[[ii]][jj,c( "low", "uppe")] <- quantile( tres[( tquants>=intev[1]) & ( tquants<intev[2])], probs = perc)
			}
			envel[[ii]] <- cbind( envel[[ii]], qnorm( envel[[ii]]))
		}
		else
			envel[[ii]] <- matrix( -99, nrow=0, ncol=8)
		colnames( envel[[ii]]) <- c( "unif", "unif.obs", "unif.low", "unif.uppe", "norm",  "norm.obs", "norm.low", "norm.uppe")
		envel[[ii]][,"norm.obs"] <- ifelse( envel[[ii]][,"unif.obs"]>=1, 50, envel[[ii]][,"norm.obs"])
		envel[[ii]][,"norm.obs"] <- ifelse( envel[[ii]][,"unif.obs"]<=0, -50, envel[[ii]][,"norm.obs"])
		envel[[ii]][,"norm.low"] <- ifelse( envel[[ii]][,"unif.low"]>=1, 50, envel[[ii]][,"norm.low"])
		envel[[ii]][,"norm.low"] <- ifelse( envel[[ii]][,"unif.low"]<=0, -50, envel[[ii]][,"norm.low"])
		envel[[ii]][,"norm.uppe"] <- ifelse( envel[[ii]][,"unif.uppe"]>=1, 50, envel[[ii]][,"norm.uppe"])
		envel[[ii]][,"norm.uppe"] <- ifelse( envel[[ii]][,"unif.uppe"]<=0, -50, envel[[ii]][,"norm.uppe"])
	}
	print( "patch residual intervals calculated")

	tcross <- array( -99, c( 3, n.cats, n.cats), list( c( "obs", "low", "uppe"), paste( "cat", 1:n.cats, sep=""), paste( "cat", 1:n.cats, sep="")))
	tcross["obs",,] <- temp[["njumps"]]
	for( ii in 1:n.cats)
		for( jj in 1:n.cats){
			tcross[ c( "low", "uppe"), ii, jj] <- quantile( cross[,ii,jj], perc)
		}
	uniPts <- (1:ncol( jumpDistRes)-0.5)/ncol( jumpDistRes)
	jumpDistResNew <- matrix( NA, nrow( jumpDistRes), ncol( jumpDistRes))
	for( ii in 1:B){
		if( any( is.na( jumpDistRes[ii,]))){
			nnna <- sum( !is.na( jumpDistRes[ii,]))
			for( jj in 1:nnna){
				intev <- (1:nnna-0.5)/nnna
				pos <- which.min( abs( uniPts-intev[jj]))
				jumpDistResNew[ii,pos] <- jumpDistRes[ii,jj]
			}
		}
		else
			jumpDistResNew[ii,] <- jumpDistRes[ii,]
		
		
	}
	jumpDistRes <- jumpDistResNew
	print( "movement residual intervals calculated")

	n.jURV <- n.cats^2 - n.cats
	very.tmp <- expand.grid( 1:n.cats, 1:n.cats)
	temp.vec.names <- double()
	for( ii in 1:( n.cats^2))
		temp.vec.names[ii] <- as.double( paste( very.tmp[ii,1], very.tmp[ii,2], sep=""))
#		temp.vec.names[ii] <- as.double( paste( very.tmp[ii,2], very.tmp[ii,1], sep=""))
#	temp.vec <- sort( c( temp$movement[ upper.tri( temp$movement, diag=F)], temp$movement[ lower.tri( temp$movement, diag=F)]))
	temp.vec <- c( temp$movement[ upper.tri( temp$movement, diag=F)], temp$movement[ lower.tri( temp$movement, diag=F)])
	temp.vec.names <- c( temp.vec.names[upper.tri( temp$movement, diag=F)], temp.vec.names[lower.tri( temp$movement, diag=F)])
	ord <- order( temp.vec)
	temp.vec <- c( (1:n.jURV-0.5)/n.jURV, temp.vec[ord], apply( jumpDistRes, MARG=2, FUN=quantile, probs=perc[1], na.rm=TRUE), apply( jumpDistRes, MARG=2, FUN=quantile, probs=perc[2], na.rm=TRUE))#, temp.vec.names[ord])
	jumpDist1 <- matrix( temp.vec, nrow=n.cats^2-n.cats, ncol=4)
	jumpDist1 <- cbind( jumpDist1, qnorm( jumpDist1))
	jumpDist1 <- cbind( jumpDist1, temp.vec.names[ord])
	colnames( jumpDist1) <- c( "unif", "unif.obs", "unif.low", "unif.uppe", "norm", "norm.obs", "norm.low", "norm.uppe", "labels")
	jumpDist1[,"norm.obs"] <- ifelse( jumpDist1[,"unif.obs"]>=1, 50, jumpDist1[,"norm.obs"])
	jumpDist1[,"norm.obs"] <- ifelse( jumpDist1[,"unif.obs"]<=0, -50, jumpDist1[,"norm.obs"])
	jumpDist1[,"norm.low"] <- ifelse( jumpDist1[,"unif.low"]>=1, 50, jumpDist1[,"norm.low"])
	jumpDist1[,"norm.low"] <- ifelse( jumpDist1[,"unif.low"]<=0, -50, jumpDist1[,"norm.low"])
	jumpDist1[,"norm.uppe"] <- ifelse( jumpDist1[,"unif.uppe"]>=1, 50, jumpDist1[,"norm.uppe"])
	jumpDist1[,"norm.uppe"] <- ifelse( jumpDist1[,"unif.uppe"]<=0, -50, jumpDist1[,"norm.uppe"])
	
	
	retu <- list()
	retu[["patch"]] <- envel
	retu[["movement"]] <- jumpDist1
	retu[["njumps"]] <- tcross
#	retu[["jumpDist"]] <- temp[["jumpDists"]]
#	retu[["jumpURV"]] <- jumpDist1
#	retu[["jumpURV1"]] <- temp[["jumpURV"]]
	
	print( "exiting")
	return( retu)
}


"examplesForDiagnostics.doc" <-
function () 
{
  print( "This function does nothing. See the examples in the documentation for code to generate the examples in the diagnostics paper")

}



"hrplot" <-
function( envel.obj, ylim=NULL, ...) {

   res <- envel.obj
	allVals <- c( res[,"norm.low"], res[,"norm.uppe"], res[,"norm.obs"])-res[,"norm"]

	if( is.null( ylim))
		limY <- c( min( allVals), max( allVals))
	else
		limY <- ylim

	if( ! "labels" %in% colnames( res)){
		plot( res[,"norm"], res[,"norm.obs"]-res[,"norm"], ylim=limY, ...)
		lines( res[,"norm"], res[,"norm.low"]-res[,"norm"], pch='.')
		lines( res[,"norm"], res[,"norm.uppe"]-res[,"norm"], pch='.')
		abline( h=0,lty=2)
	}
	else {
		plot( res[,"norm"], res[,"norm.obs"]-res[,"norm"], ylim=limY, ...,type='n')
		text( res[,"norm"], res[,"norm.obs"]-res[,"norm"], labels=res[,"labels"], ...)
		lines( res[,"norm"], res[,"norm.low"]-res[,"norm"], pch='.')
		lines( res[,"norm"], res[,"norm.uppe"]-res[,"norm"], pch='.')
		abline( h=0,lty=2)
	}


}



"MVfill" <-
function ( fm, states=NULL, chain.id=NULL, X=NULL) 
{
  if( is.null( states)){
    stop( error="No observed state data supplied")
    res <- NULL
  }
  else{
    #get states into an understandable order
    states <- factor( states)
    levels( states) <- order( levels( states))
    n.cats <- length( levels( states))
    states <- matrix( as.numeric( states), ncol=1)
    #clean up chain.id input (if required)
    if( is.null( chain.id)){
      chain.id <- rep( 1, nrow( states))
      cat( "No chain.id supplied - all observations assumed to be from a single chain", "\n")
    }
    chain.id <- factor( chain.id)
    n.obs <- summary( chain.id) 
    chain.id <- matrix( as.numeric( chain.id), ncol=1)
    #clean up X matrix ( if required)   
    if( is.null( X)){   
      X <- matrix( 1, nrow=nrow( states), ncol=1)   
      cat( "No design matrix supplied - assumed only to be the constant", "\n")
    }
    else
      X <- data.matrix( X)
    #check if fm is given
    if( is.null( fm))
    	stop( "No model given - function aborted", "\n")
    #check if states, chains and X agree
    if( ( nrow( states) != nrow( chain.id)) || ( nrow( states) != nrow( X)) || ( nrow( chain.id) != nrow( X)))
      stop( "Error: states, chain.id, and/or X have different lengths")
    pars <- fm$pars
    phi.id <- fm$stuff$phi.id
    pi.id <- fm$stuff$pi.id
    if( ( nrow( phi.id) != ncol( X)) || ( length( pi.id) != ncol( X)))
    	stop( "X and fm do not agree", "\n")
    #top and tailing each chain for NAs
    id.vec <- rep( 0, nrow( states))
    states <- ifelse( is.na( states), -9, states)
    id.vec <- .C( "topTailChains", as.integer( states), as.integer( n.obs), as.integer( length( n.obs)), idvec=as.integer( id.vec),PACKAGE="RMC")$idvec
    states.red <- states[id.vec!=1,,drop=FALSE]
    chain.id.red <- chain.id[id.vec!=1,,drop=FALSE]
    X.red <- X[id.vec!=1,,drop=FALSE]
    n.obs <- summary( as.factor( chain.id.red))
    #Defining and filling in the conditional probabilities
    temp <- .C( "MVfillCpp", temp=double( length( states.red)*n.cats), as.double( pars), as.integer(n.cats), 
    		as.integer( ncol( X.red)), as.integer( nrow( X.red)), as.integer( length( n.obs)), as.double( X.red), 
    		as.integer( phi.id), as.integer( n.obs), as.integer( states.red), as.integer( pi.id),PACKAGE="RMC")$temp
    probs <- matrix( NA, nrow=nrow( X), ncol=n.cats)
    probs[id.vec!=1,] <- temp
  }
  return( probs)
}



"nd2" <-
function(x0, f, m=NULL, D.accur=4, eps=NULL, ...) {
# A function to compute highly accurate first-order derivatives 
# Stolen (mostly) from the net and adapted / modified by Scott
# From Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)  

# x0 is the point where the derivative is to be evaluated, 
# f is the function that requires differentiating
# m is output dimension of f, that is f:R^n -> R^m
#D.accur is the required accuracy of the resulting derivative. Options are 2 and 4. The 2 choice does a two point finite difference approximation and the 4 choice does a four point finite difference approximation.
#eps is the 

# Report any bugs to Scott as he uses this extensively!
  D.n<-length(x0)
  if (is.null(m)) {
    D.f0<-f(x0, ...)
    m<-length(D.f0) }
  if (D.accur==2) {
    D.w<-tcrossprod(rep(1,m),c(-1/2,1/2))
    D.co<-c(-1,1) }
  else {
    D.w<-tcrossprod(rep(1,m),c(1/12,-2/3,2/3,-1/12))
    D.co<-c(-2,-1,1,2) }
  D.n.c<-length(D.co)
  if( is.null( eps)) {
    macheps<-.Machine$double.eps
    D.h<-macheps^(1/3)*abs(x0)
  }
  else
    D.h <- rep( eps, D.accur)
  D.deriv<-matrix(NA,nrow=m,ncol=D.n)
  for (ii in 1:D.n) {
    D.temp.f<-matrix(0,m,D.n.c)
    for (jj in 1:D.n.c) {
      D.xd<-x0+D.h[ii]*D.co[jj]*(1:D.n==ii)
      D.temp.f[,jj]<-f(D.xd, ...) }
    D.deriv[,ii]<-rowSums(D.w*D.temp.f)/D.h[ii] }
  return( D.deriv)
}


".onLoad" <-
function( libname, pkgname){
   # Generic DLL loader
   dll.path <- file.path( libname, pkgname, 'libs')
   if( nzchar( subarch <- .Platform$r_arch))
     dll.path <- file.path( dll.path, subarch)
   this.ext <- paste( sub( '.', '[.]', .Platform$dynlib.ext, fixed=TRUE), '$', sep='')

   dlls <- dir( dll.path, pattern=this.ext, full=FALSE)
   names( dlls) <- dlls
   if( length( dlls)) 
     lapply( dlls, function( x) library.dynam( sub( this.ext, '', x), package=pkgname, lib.loc=libname))
}


"RMC.GoF" <-
function ( fm, states=NULL, chain.id=NULL, X=NULL, lambda=2/3) 
{
  if( is.null( states)){
    stop( error="No observed state data supplied")
    res <- NULL
  }
  else{
    #get states into an understandable order
    states <- factor( states)
    levels( states) <- order( levels( states))
    n.cats <- length( levels( states))
    states <- matrix( as.numeric( states), ncol=1)
    #clean up chain.id input (if required)
    if( is.null( chain.id)){
      chain.id <- rep( 1, nrow( states))
      cat( "No chain.id supplied - all observations assumed to be from a single chain", "\n")
    }
    chain.id <- factor( chain.id)
    n.obs <- summary( chain.id) 
    chain.id <- matrix( as.numeric( chain.id), ncol=1)
    #clean up X matrix ( if required)   
    if( is.null( X)){   
      X <- matrix( 1, nrow=nrow( states), ncol=1)   
      cat( "No design matrix supplied - assumed only to be the constant", "\n")
    }
    else
      X <- data.matrix( X)
    #check if fm is given
    if( is.null( fm))
    	stop( "No model given - function aborted", "\n")
    #check if states, chains and X agree
    if( ( nrow( states) != nrow( chain.id)) || ( nrow( states) != nrow( X)) || ( nrow( chain.id) != nrow( X)))
      stop( "Error: states, chain.id, and/or X have different lengths")
    pars <- fm$pars
    phi.id <- fm$stuff$phi.id
    pi.id <- fm$stuff$pi.id
    if( ( nrow( phi.id) != ncol( X)) || ( length( pi.id) != ncol( X)))
    	stop( "X and fm do not agree", "\n")
    #top and tailing each chain for NAs
    id.vec <- rep( 0, nrow( states))
    states <- ifelse( is.na( states), -9, states)
    id.vec <- .C( "topTailChains", as.integer( states), as.integer( n.obs), as.integer( length( n.obs)), idvec=as.integer( id.vec),PACKAGE="RMC")$idvec
    states.red <- states[id.vec!=1,,drop=FALSE]
    chain.id.red <- chain.id[id.vec!=1,,drop=FALSE]
    X.red <- X[id.vec!=1,,drop=FALSE]
    n.obs <- summary( as.factor( chain.id.red))
    #Defining and filling in the conditional probabilities
    temp <- .C( "getCondPreds", temp=double( length( states.red)*n.cats), as.double( pars), as.integer(n.cats), 
    		as.integer( ncol( X.red)), as.integer( nrow( X.red)), as.integer( length( n.obs)), as.double( X.red), 
    		as.integer( phi.id), as.integer( n.obs), as.integer( states.red), as.integer( pi.id),PACKAGE="RMC")$temp
    probs <- matrix( NA, nrow=nrow( X), ncol=n.cats)
    probs[id.vec!=1] <- temp
    probs[probs==-9] <- NA

    ustat.vec <- calc.ustat.vec(states, probs, lambda)
    estat.vec <- calc.estat.vec(probs, lambda)
    mart <- ustat.vec - estat.vec
    In.stat <- ( 2/( lambda*(lambda+1))) * sum( mart, na.rm=T)

#    xistat.vec <- calc.xistat.vec( states, probs, lambda, fm$vcov)

  }

  return( NULL)

}


"RMC.mod" <-
function( states, chain.id=NULL, X=NULL, phi.id=NULL, pi.id=NULL, vcov=FALSE, inits=NULL, contr=list( maxit=1000, epsg=1e-8, epsf=1e-8, epsx=1e-8), quiet=FALSE, penalty=0)
{
#All ploglike calculated using C function with a few speed tricks
#penalty argument is the variance of a normal penalty on the parameters.  Only use(?) when numerics are an issue.

  if( is.null( states)){
    stop( error="No observed state data supplied")
    res <- NULL
  }
  else{
    #get states into an understandable order
    states <- factor( states)
    levels( states) <- as.character( order( as.numeric( levels( states))))
    n.cats <- length( levels( states))
    states <- matrix( as.numeric( states), ncol=1)
    #clean up chain.id input (if required)
    if( is.null( chain.id)){
      chain.id <- rep( 1, nrow( states))
      if( !quiet)
        cat( "No chain.id supplied - all observations assumed to be from a single chain", "\n")
    }
    chain.id <- factor( chain.id)
    n.obs <- summary( chain.id) 
    chain.id <- matrix( as.numeric( chain.id), ncol=1)
    #clean up X matrix ( if required)   
    if( is.null( X)){   
      X <- matrix( 1, nrow=nrow( states), ncol=1)   
      if( !quiet)
        cat( "No design matrix supplied - assumed only to be the constant", "\n")
    }
    else
      X <- data.matrix( X)
    #clean up phi.id and pi.id if required
    if( is.null( phi.id) || nrow( phi.id) != ncol( X)){
      phi.id <- matrix( 1, nrow=ncol( X), ncol=n.cats)
      if( !quiet)
        cat( "Incidence matrix for phi not given or doesn't match design matrix - all variables included", "\n")
    }
    if( is.null( pi.id) || length( pi.id) != ncol( X)){
      pi.id <- rep( 1, ncol( X))
      if( !quiet) 
        cat( "Incidence vector for pi not given or doesn't match design matrix - all variables included", "\n")
    }
    #top and tailing each chain for NAs
    id.vec <- rep( 0, nrow( states))
    states <- ifelse( is.na( states), -9, states)
    id.vec <- .C( "topTailChains", as.integer( states), as.integer( n.obs), as.integer( length( n.obs)), idvec=as.integer( id.vec),PACKAGE="RMC")$idvec
    states <- states[id.vec!=1,,drop=FALSE]
    #check to see if there are observed states with NA covariates -- presumably because of filled covariates from MVfill
    tmp.X.NA <- is.na( rowSums( X))#apply( X, 1, function(x) any( is.na( x)))
    tmp.which <- which( tmp.X.NA[id.vec!=1] & !is.na( states))
    if( length( tmp.which != 0)){
      print( "Exitting prior to estimation -- there are missing values in the covariates where they ought not to be")
      return()
    }
    chain.id <- chain.id[id.vec!=1,,drop=FALSE]
    X <- X[id.vec!=1,,drop=FALSE]
    n.obs <- summary( as.factor( chain.id))
    #preparing parameter vector etc.
    if( is.null( inits) || length( inits)!= ( sum( phi.id) + sum( pi.id)*( n.cats-1))){
      pars <- rep( 0, sum( phi.id) + sum( pi.id)*( n.cats-1))
      if( !quiet)
        cat( "All initial values assigned to zero", "\n")
    }
    else
      pars <-inits
      
    #Doing the estimation
    if( !quiet){
      cat( "\n", "Performing Estimation of Parameters for Reversible Markov Model", "\n")
      cat( "\n", "Convergence trace", "\n")
      cat( "  iter  ", "-pllike: ", paste( "phi", 1:sum( phi.id), sep="_", collapse="    "), "   ")
      for( ii in 1:sum( pi.id))
        cat( paste( paste( "pi", ii, sep="_"), 1:( n.cats-1), sep="", collapse="    "), "   ")
      cat("\n")
    }
    #Fixing up control parameters for L-BFGS
    contr$mNo <- contr$maxit
                                                                                                                                                                                
    RMC.time <- double( 2)
    RMC.time <- 0
    names( RMC.time) <- "estimation"
    RMC.time["estimation"] <- system.time( 
    	temp <- .C( "RMC", pars=as.double( pars), plogl=double( 1), scores=double( length( pars)), hessian=double( length( pars)^2), as.integer( vcov),
    		conv=integer( 1), iters=integer( 1), as.integer( n.cats), as.integer( length( pi.id)), as.integer( length( states)), 
    		as.integer( length( n.obs)), as.double( X), as.integer( phi.id), as.integer( n.obs), as.integer( states), as.integer( pi.id), 
    		as.integer( contr$mNo), as.double( contr$epsg), as.double( contr$epsf), as.double( contr$epsx), as.integer( contr$maxit), as.integer( quiet), as.double( penalty),PACKAGE="RMC"))[1]

   #COV stuff                             
   if( vcov){
	if( !quiet)
          cat( "\n", "Calculating Variances for the Parameter Estimates", "\n\n")
    	vcov.mat <- try( solve( - matrix( temp$hessian, ncol=length( pars))))
    	if( inherits( vcov.mat, "try-error", FALSE))
        vcov.mat <- matrix( 0*temp$hessian, ncol=length( pars))
   }
	else
		vcov.mat <- matrix( 1, ncol=1, nrow=1)
	if( !quiet){
          cat( " Computing Done ")
          if( temp$conv==-1)
            cat( "-- Error in code")
          if( temp$conv==5)
            cat( "-- not converged")
          cat( "\n")
        }

	stuff <- list( phi.id=phi.id, pi.id=pi.id, n.cats=n.cats, n.covars=ncol( X), n.chains=length( n.obs), n.obs=n.obs, n.pars=length( temp$pars), id.vec=id.vec)

	res <- list( pars=temp$pars, like=temp$plogl, scores=temp$scores, vcov=vcov.mat, conv=temp$conv, time=RMC.time, niter=temp$iter, stuff=stuff)
	}                                                                                                                                                                             
	return( res)
}



"RMC.pred" <-
function ( fit, fit2=NULL, pts) 
{
	if( length( as.double( pts)) == length( fit$stuff$pi.id))
	  pts <- matrix( pts, nrow=1)
   if( is.data.frame( pts))
     pts <- as.matrix( pts)
   if( length( as.double( pts)) == length( fit$stuff$pi.id))
     pts <- matrix( pts, nrow=1)
              		
	####	Single level of prediction
	if( is.null( fit2)){

     if( all( fit$vcov!=0 && fit$vcov!=1))
       do.cov <- 1
     else
       do.cov <- 0
  
     temp <- .C( "predPtRMC", Rpars=as.double( fit$pars), Rpts=as.double( pts), Rncats=as.integer( fit$stuff$n.cats), Rncovars=as.integer( fit$stuff$n.covars), Rntot=as.integer( nrow( pts)), RphiID=as.integer( fit$stuff$phi.id), RpiID=as.integer( fit$stuff$pi.id), RdumbStates=as.integer( rep( 1, nrow( pts))), ppts=double( nrow( pts)*fit$stuff$n.cats), area=double( fit$stuff$n.cats), Rjac=double( length( fit$pars) * fit$stuff$n.cats), RdoCov=as.integer( do.cov),PACKAGE="RMC")
          
     ppts <- matrix( temp$ppts, nrow=nrow(pts), ncol=fit$stuff$n.cats)
     area.pred <- temp$area
  
     if( do.cov==1){
       jac <- matrix( temp$Rjac, nrow=fit$stuff$n.cats, ncol=length( fit$pars))  
       var <- jac%*%tcrossprod(fit$vcov,jac)
     }
     else
     var=NULL
  
     return( list( area=area.pred, vcov=var, pts=ppts))
  }  
  ####	2 levels of prediction
  else{
		pptsF <- matrix( -99, nrow=nrow( pts)*fit$stuff$n.cats, ncol=fit$stuff$n.cats)
		pptsG <- matrix( -99, nrow=nrow( pts), ncol=fit2$stuff$n.cats)
		FptsPreds <- matrix( -99, nrow=nrow(pts), ncol=fit$stuff$n.cats)
		
		indMat <- cbind( rep( 1, fit2$stuff$n.cats), diag( fit2$stuff$n.cats))
		ptsF <- kronecker( indMat, pts)
		
		area.fun <- function(tmp.pars, do.pts=FALSE){
 			do.cov.tmp <- 0
			G.pars <- tmp.pars[1:length(fit2$pars)]
			F.pars <- tmp.pars[length(fit2$pars)+1:length(fit$pars)]
			tempG <- .C( "predPtRMC", Rpars=as.double( G.pars), Rpts=as.double( pts), Rncats=as.integer( fit2$stuff$n.cats), Rncovars=as.integer( fit2$stuff$n.covars), Rntot=as.integer( nrow( pts)), RphiID=as.integer( fit2$stuff$phi.id), RpiID=as.integer( fit2$stuff$pi.id), RdumbStates=as.integer( rep( 1, nrow( pts))), ppts=double( nrow( pts)*fit2$stuff$n.cats), area=double( fit2$stuff$n.cats), Rjac=double( length( fit2$pars) * fit2$stuff$n.cats), RdoCov=as.integer( do.cov.tmp),PACKAGE="RMC")
			tempF <- .C( "predPtRMC", Rpars=as.double( F.pars), Rpts=as.double( ptsF), Rncats=as.integer( fit$stuff$n.cats), Rncovars=as.integer( fit$stuff$n.covars), Rntot=as.integer( nrow( ptsF)), RphiID=as.integer( fit$stuff$phi.id), RpiID=as.integer( fit$stuff$pi.id), RdumbStates=as.integer( rep( 1, nrow( ptsF))), ppts=double( nrow( ptsF)*fit$stuff$n.cats), area=double( fit$stuff$n.cats), Rjac=double( length( fit$pars) * fit$stuff$n.cats), RdoCov=as.integer( do.cov.tmp),PACKAGE="RMC")
			pptsG <- matrix( tempG$ppts, nrow=nrow(pts), ncol=fit2$stuff$n.cats)
			pptsF <- matrix( tempF$ppts, nrow=nrow(ptsF), ncol=fit$stuff$n.cats)
			pptsF <- pptsF * as.double( pptsG)

			area.pred1 <- colSums( pptsF) / nrow( pts)							
			if( do.pts){
				tmp <- rep( 1:nrow( pts), fit2$stuff$n.cats)
				FptsPreds <- matrix( -99, nrow=nrow( pts), ncol=fit$stuff$n.cats)
				for( ii in 1:fit$stuff$n.cats)
					FptsPreds[,ii] <- tapply( pptsF[,ii], tmp, sum)
				return( list( pts=FptsPreds, area=area.pred1))
			}
			else
				return( area.pred1)
		}			

		if( all( fit$vcov!=0 && fit$vcov!=1 && fit2$vcov!=0 && fit2$vcov!=1))
			do.cov <- 1
		else
	      do.cov <- 0

		if( do.cov==1){
			jac <- nd2( x0=c(fit2$pars, fit$pars), f=area.fun)
			ng <- length( fit2$pars)
			nf <- length( fit$pars)
			tot.var <- matrix( 0, nrow=ng+nf, ncol=ng+nf)
			tot.var[1:ng, 1:ng] <- fit2$vcov
			tot.var[ng+1:nf, ng+1:nf] <- fit$vcov
		
			vari <- jac %*% tcrossprod( tot.var, jac)
		}
		else
			vari <- NULL

		area.pred <- area.fun( c(fit2$pars, fit$pars), do.pts=TRUE)

		return( list( area=area.pred$area, vcov=vari, pts=area.pred$pts))
	}	
}



"sim.chain" <-
function (n.chains = 5, n.obs = rep( 100, n.chains), n.cats = 3, n.covars = 1, beta = NULL, gamma = NULL, X = NULL) 
{
   if (n.covars < 1) {
      n.covars <- 1
      print( "inappropriate number of chains, changing to n.chains=1")
   }
   dat <- matrix( -99, ncol = 2 + n.covars, nrow = sum(n.obs))
   if (n.covars == 1) 
      colnames(dat) <- c("chain", "state", "const")
   else {
      if (is.null(X)) 
         colnames(dat) <- c("chain", "state", "const", paste("var", 1:(n.covars - 1), sep = ""))
     else colnames(dat) <- c("chain", "state", colnames(X))
   }
   place <- 0
   for (ii in 1:n.chains) {
      dat[place + 1:n.obs[ii], "chain"] <- ii
      place <- place + n.obs[ii]
   }
   if (is.null(X)) {
      dat[, "const"] <- rep(1, sum(n.obs))
      if (n.covars != 1) 
         dat[, -(1:3)] <- scale(rnorm(n = sum(n.obs) * (n.covars - 1)))
   }
   else{
   	dat[, -(1:2)] <- data.matrix( X)
   }
   if (is.null(beta) | is.null(gamma)) {
      beta <- matrix(rnorm(n = n.covars * n.cats), ncol = n.cats, nrow = n.covars)
      beta[, 1] <- 0
      gamma <- matrix(rnorm(n = n.covars * n.cats), ncol = n.cats, nrow = n.covars)
      if (n.covars > 1) 
         rownames(beta) <- rownames(gamma) <- c("const", paste("covar", 1:(n.covars - 1), sep = ""))
      colnames(beta) <- paste("categ", 1:n.cats, sep = "")
      colnames(gamma) <- paste("categ", 1:n.cats, sep = "")
   }
   else
   	if( is.double( beta))
   		beta <- matrix( beta, nrow=n.covars, ncol=n.cats)
   
   pi.id <- rep( 1, n.covars)
   phi.id <- matrix( 1, nrow=n.covars, ncol=n.cats)
   paras <- c( as.double( gamma), as.double( beta[,-1]))
   
	temp <- .C( "simChains",	states=as.integer( dat[,2]), as.integer( n.obs), as.double( dat[,-(1:2)]), as.integer( n.cats), as.integer( n.covars), as.integer( sum( n.obs)), as.integer( length( n.obs)), as.double( paras), as.integer( phi.id), as.integer( pi.id), marg=as.double( rep( 0, n.cats)),PACKAGE="RMC")
	
	dat[,"state"] <- temp$states
#	attr(dat, "statDist") <- temp$marg
	   
	return( dat)
}



"simRandWalk" <-
function ( nc=5, ni=rep( 1000, nc), init.var=1, seq.var=0.1) 
{
	res <- matrix( -99, ncol=2, nrow=sum( ni))
	colnames( res) <- c( "chain", "var1")
	place <- 1
	for( ii in 1:nc){
		res[place, "chain"] <- ii
		res[place, "var1"] <- rnorm( 1, sd=sqrt( init.var)		)
		place <- place+1
		for( jj in 2:ni[ii]){
			res[place, "chain"] <- ii
			res[place, "var1"] <- res[place-1, "var1"] + rnorm( 1, sd=sqrt( seq.var))
			place <- place + 1
		}
	}	
	return( res)
}


