# This is package fishMod 

"deltaLN" <-
function ( ln.form, binary.form, data, residuals=TRUE) 
{
  temp <- model.frame( ln.form, data=as.data.frame( data))
  X <- model.matrix( ln.form, temp)
  offy <- model.offset( temp)
  if( is.null( offy))
    offy <- rep( 0, nrow( X))

  nonzeros <- model.response( temp)>0
  temp.ln <- model.frame( ln.form, data=as.data.frame( data[nonzeros,]))
  nz.y <- log( model.response( temp.ln))
  nz.X <- model.matrix( ln.form, temp.ln)
  nz.offset <- model.offset( temp.ln)
  if (is.null(nz.offset)) 
    nz.offset <- rep(0, length(nz.y))

  temp.bin <- model.frame( binary.form, data=as.data.frame( data))
  temp.bin <- model.matrix( binary.form, temp.bin)
  bin.X <- cbind( as.numeric( nonzeros), temp.bin)
  colnames( bin.X) <- c("positive", colnames( bin.X)[-1])
#  if( length( colnames( temp.bin))==1)
#    binary.form <- as.formula( paste( colnames( temp.bin)[1],"~1"))
#  else
#    binary.form <- as.formula( paste( colnames( temp.bin)[1],"~1+",paste( colnames( temp.bin)[-1], collapse="+")))

  #positive data log-normal
  lnMod <- lm( nz.y~-1+nz.X, offset=nz.offset)#ln.form, temp.ln)
  #bianry glm
  binMod <- glm( bin.X[,"positive"]~-1+bin.X[,-1], family=binomial())

  stdev <- summary( lnMod)$sigma

  coefs <- list( binary=binMod$coef, ln=lnMod$coef, ln.sigma2=stdev^2)

  logl <- sum( dlnorm( exp(nz.y), lnMod$fitted, sdlog=stdev, log=TRUE)) + sum( dbinom( nonzeros, size=1, prob=binMod$fitted, log=TRUE))
  n <- nrow( temp)
  ncovars <- length( lnMod$coef) + 1 + length( binMod$coef)
  nnonzero <- sum( nonzeros)
  nzero <- n-nnonzero
  AIC <- -2*logl + 2*ncovars
  BIC <- -2*logl + log( n)*ncovars

  fitted <- var <- rep( -999, n)
#  colnames( temp)[ncol( temp)] <- "nz.offset"
#  lpv <- predict( lnMod, temp)
  lpv <- X %*% lnMod$coef + offy
  pv <- exp( lpv + 0.5*(stdev^2))
  fitted <- binMod$fitted * pv
	var <- binMod$fitted*exp( 2*lpv+stdev^2)*(exp( stdev^2)-binMod$fitted)
  
  if( residuals){
    resids <- matrix( rep( -999, 2*n), ncol=2, dimnames=list(NULL,c("quantile","Pearson")))
    resids[nonzeros,"quantile"] <- (1-binMod$fitted[nonzeros]) + binMod$fitted[nonzeros] * plnorm( temp.ln[,1], lnMod$fitted, stdev, log.p=FALSE)
    resids[!nonzeros,"quantile"] <- runif( n=sum( !nonzeros), min=rep( 0, nzero), max=1-binMod$fitted[!nonzeros])
    resids[,"quantile"] <- qnorm( resids[,"quantile"])
    resids[,"Pearson"] <- ( model.response(temp)-fitted) / sqrt(var)
  }
  else
    resids <- NULL
  

  res <- list( coefs=coefs, logl=logl, AIC=AIC, BIC=BIC, fitted=fitted, fittedVar=var, residuals=resids, n=n, ncovars=ncovars, nzero=nzero, lnMod=lnMod, binMod=binMod)

  class( res) <- "DeltaLNmod"

  return( res)
}


"dPoisGam" <-
function ( y, lambda, mu.Z, alpha, LOG=TRUE) 
{
#function to calculate Random sum (Tweedie) densities.
#y is the value of the r.v.  Can be a vector
#mu.N is the mean of the Poisson summing r.v. Can be a vector of length(y)
#mu.Z is the mean of the Gamma rv Can be a vector of length(y)
#alpha is the `other' parameter of the gamma distribution s.t. var = ( mu.Z^2)/alpha Can be a vector of length(y)
#If mu.N, mu.Z or alpha are scalare but y isn't then they will be used for all y. If lengths mis-match then error
#LOG=TRUE gives the density on the log scale
#do.checks=TRUE checks the input vectors for compatability and gives errors / changes them as appropriate.
#do.checks=FALSE doesn't check and relies on the user to have things right. If not right then catastrophic failure may occur.

#  if( any( is.null( c( y, mu.N, mu.Z, alpha)))){
#    print( "Error: null input values -- please check.  Null values are:")
#    tmp <- double( is.null( c( y, mu.N, mu.Z, alpha)))
#    names( tmp) <- c( "y", "mu.N","mu.Z","alpha")
#    print( tmp)
#    print( "Exitting")
#    return()
#  }
  mu.N <- lambda
  if( !all( is.element( c( length( mu.N), length( mu.Z), length( alpha)), c( length( y), 1)))){
    print( "Error: length of parameter vectors does not match length of random variable vector")
    print( "Exitting")
    return()
  }

  if( length( mu.N) != length( y))
    mu.N <- rep( mu.N, length( y))
  if( length( mu.Z) != length( y))
    mu.Z <- rep( mu.Z, length( y))
  if( length( alpha) != length( y))
    alpha <- rep( alpha, length( y))

  res <- .Call( "dTweedie", as.numeric( y), as.numeric( mu.N), as.numeric( mu.Z), as.numeric( alpha), as.integer( LOG), PACKAGE="fishMod")

  return( res)

}


"dPoisGamDerivs" <-
function ( y=NULL, lambda=NULL, mu.Z=NULL, alpha=NULL, do.checks=TRUE) 
{
#function to calculate Random sum (Tweedie) densities.
#y is the value of the r.v.  Can be a vector
#mu.N is the mean of the Poisson summing r.v. Can be a vector of length(y)
#mu.Z is the mean of the Gamma rv Can be a vector of length(y)
#alpha is the `other' parameter of the gamma distribution s.t. var = ( mu.Z^2)/alpha Can be a vector of length(y)
#If mu.N, mu.Z or alpha are scalare but y isn't then they will be used for all y. If lengths mis-match then error
#LOG=TRUE gives the density on the log scale
#do.checks=TRUE checks the input vectors for compatability and gives errors / changes them as appropriate.
#do.checks=FALSE doesn't check and relies on the user to have things right. If not right then catastrophic failure may occur.

  mu.N <- lambda
  if( do.checks){
    if( any( is.null( c( y, mu.N, mu.Z, alpha)))){
      print( "Error: null input values -- please check.  Null values are:")
      tmp <- double( is.null( c( y, mu.N, mu.Z, alpha)))
      names( tmp) <- c( "y", "mu.N","mu.Z","alpha")
      print( tmp)
      print( "Exitting")
      return()
    }

    if( !all( is.element( c( length( mu.N), length( mu.Z), length( alpha)), c( length( y), 1)))){
      print( "Error: length of parameter vectors does not match length of random variable vector")
      print( "Exitting")
    }

    if( length( mu.N) != length( y))
      mu.N <- rep( mu.N, length( y))
    if( length( mu.Z) != length( y))
      mu.Z <- rep( mu.Z, length( y))
    if( length( alpha) != length( y))
      alpha <- rep( alpha, length( y))
  }

  res <- .Call( "dTweedieDeriv", as.numeric( y), as.numeric( mu.N), as.numeric( mu.Z), as.numeric( alpha), PACKAGE="fishMod")
  colnames( res) <- c("lambda","mu.Z","alpha")
  return( res)

}


"dTweedie" <-
function ( y, mu, phi, p, LOG=TRUE) 
{
  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*mu^(p-1)
  mu.Z <- alpha * tau

  dens <- dPoisGam( y, lambda, mu.Z, alpha, LOG)
  
  return( dens)

}


"ldPoisGam.lp" <-
function ( parms, y, X.p, X.g, offsetty, alpha, wts=rep( 1, length( y))) 
{
  mu.p <- exp( X.p %*% parms[1:ncol( X.p)] + offsetty)
  mu.g <- exp( X.g %*% parms[ncol( X.p)+1:ncol( X.g)])

  if( is.null( alpha))
    return( - sum( wts * dPoisGam( y, lambda=mu.p, mu.Z=mu.g, alpha=exp(tail( parms, 1)), LOG=TRUE)))
  else
    return( - sum( wts * dPoisGam( y, lambda=mu.p, mu.Z=mu.g, alpha=alpha, LOG=TRUE)))
}


"ldPoisGam.lp.deriv" <-
function ( parms, y, X.p, X.g, offsetty, alpha, wts=rep( 1, length( y))) 
{
  mu.p <- exp( X.p %*% parms[1:ncol( X.p)] + offsetty)
  mu.g <- exp( X.g %*% parms[ncol( X.p)+1:ncol( X.g)])

  if( is.null( alpha))
    alpha1 <- exp( tail( parms,1))
  else
    alpha1 <- alpha

  dTweedparms <- - wts * dPoisGamDerivs( y, lambda=mu.p, mu.Z=mu.g, alpha=alpha1)		#wts should replicate appropriately

  deri.lambda <- ( as.numeric( dTweedparms[,"lambda"] * mu.p)) * X.p
  deri.mu <- ( as.numeric( dTweedparms[,"mu.Z"] * mu.g)) * X.g
  deri.alpha <- as.numeric( dTweedparms[,"alpha"]) * alpha1

  deri.all <- c( colSums( deri.lambda), colSums( deri.mu), sum( deri.alpha))

  if( is.null( alpha))
    return( deri.all)
  else
    return( deri.all[-length( deri.all)])

}


"ldTweedie.lp" <-
function ( parms, y, X.p, offsetty, phi, p, wts=rep( 1, length( y))) 
{
  mu <- exp( X.p %*% parms[1:ncol( X.p)] + offsetty)

  if( is.null( phi) & is.null( p)){
    phi <- parms[ncol( X.p) + 1]
    p <- parms[ncol( X.p) + 2]
  }
  if( is.null( phi) & !is.null( p))
    phi <- parms[ncol( X.p)+1]
  if( !is.null( phi) & is.null( p))
    p <- parms[ncol( X.p)+1]
  
  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*(mu^(p-1))
  mu.Z <- alpha * tau

  return( -sum( wts * dPoisGam( y, lambda=lambda, mu.Z=mu.Z, alpha=alpha, LOG=TRUE)))
}


"ldTweedie.lp.deriv" <-
function ( parms, y, X.p, offsetty, phi, p, wts=rep( 1, length( y))) 
{
  mu <- exp( X.p %*% parms[1:ncol( X.p)] + offsetty)

  p.flag <- phi.flag <- FALSE
  if( is.null( phi) & is.null( p)){
    p.flag <- phi.flag <- TRUE
    phi <- parms[ncol( X.p) + 1]
    p <- parms[ncol( X.p) + 2]
  }
  if( is.null( phi) & !is.null( p)){
    phi <- parms[ncol( X.p)+1]
    phi.flag <- TRUE
  }
  if( !is.null( phi) & is.null( p)){
    p <- parms[ncol( X.p)+1]
    p.flag <- TRUE
  }

  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*mu^(p-1)
  mu.Z <- alpha * tau

  dTweedparms <- -wts * dPoisGamDerivs( y, lambda=lambda, mu.Z=mu.Z, alpha=alpha)

  DTweedparmsDmu <- matrix( c( ( mu^(1-p)) / phi, alpha*phi*( ( p-1)^2)*( mu^(p-2)), rep( 0, length( mu))), nrow=3, byrow=T)
  tmp <- rowSums( dTweedparms * t( DTweedparmsDmu))
  tmp <- tmp * mu
  tmp <- apply( X.p, 2, function( x) x*tmp)
  
  derivs <- colSums( tmp)
  
  if( phi.flag){
    DTweedparmsDphi <- matrix( c( -( ( mu^(2-p)) / ( ( phi^2)*(2-p))), alpha*( p-1)*( mu^( p-1)), rep( 0, length( mu))), nrow=3, byrow=T)
    tmpPhi <- rowSums( dTweedparms * t( DTweedparmsDphi))	#vectorised way of doing odd calculation
    derivs <- c( derivs, sum( tmpPhi))
    names( derivs)[length( derivs)] <- "phi"
  }
  if( p.flag){
    dalphadp <- -( 1+alpha) / ( p-1)
    DTweedparmsDp <- matrix( c( lambda*( 1/(2-p) - log( mu)), mu.Z*( dalphadp/alpha + 1/( p-1) + log( mu)), rep( dalphadp, length( y))), nrow=3, byrow=T)
    tmpP <- rowSums( dTweedparms * t( DTweedparmsDp))
    derivs <- c( derivs, sum( tmpP))
    names( derivs)[length( derivs)] <- "p"
  }
  
  return( derivs)
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

   dlls <- dir( dll.path, pattern=this.ext, full.names=FALSE)
   names( dlls) <- dlls
   if( length( dlls)) 
     lapply( dlls, function( x) library.dynam( sub( this.ext, '', x), package=pkgname, lib.loc=libname))
}


"pgm" <-
function ( p.form, g.form, data, wts=NULL, alpha=NULL, inits=NULL, vcov=TRUE, residuals=TRUE, trace=1)
{
  if( is.null( wts))
    wts <- rep( 1, nrow( data))
#  temp.p <- model.frame( p.form, data=as.data.frame( data))#, weights=wts)
  e<-new.env()
  e$wts <- wts
  environment(p.form) <- e
  temp.p <- model.frame( p.form, data=as.data.frame( data), weights=wts)

  y <- model.response( temp.p)
  names( y) <- NULL
  X.p <- model.matrix( p.form, data)
  offset.p <- model.offset( temp.p)
  if ( is.null( offset.p))
    offset.p <- rep( 0, nrow( temp.p))
  wts1 <- as.vector( model.weights( temp.p))
#  if( is.null( wts))
#    wts <- rep( 1, length( y))

  tmp.form <- g.form
  tmp.form[[2]] <- p.form[[2]]
  tmp.form[[3]] <- g.form[[2]]
  temp.g <- model.frame( tmp.form, as.data.frame( data))
  names( y) <- NULL
  X.g <- model.matrix( tmp.form, data)

#  if( is.null( weights)){
#    weights <- rep( 1, length( y))
#    if( trace!=0)
#      print( "all weights are unity")
#  }
#  if( length( weights) != length( y)){
#    print( "number of weights does not match number of observations")
#    return( NULL)
#  }
  
  if( is.null( inits) & is.null( alpha))
    inits <- rep( 0, ncol( X.p)+ncol( X.g)+1)
  if( is.null( inits) & !is.null( alpha))
    inits <- rep( 0, ncol( X.p)+ncol( X.g))
  if( !is.element( length( inits),ncol( X.p)+ncol( X.g)+0:1)) {
    print( "Initial values supplied are of the wrong length -- please check")
    tmp <- c( length( inits), ncol( X.p)+ncol( X.g)+1)
    names( tmp) <- c( "inits", "ncolDesign")
    return( tmp)
  }

  if( trace!=0){
    print( "Estimating parameters")
    if( is.null( alpha))
      cat("iter:", "-logl", paste( colnames(X.p),"Poisson",sep='.'), paste( colnames(X.g),"Gamma",sep='.'), "log( alpha)", "\n", sep = "\t") 
    else
      cat("iter:", "-logl", paste( colnames(X.p),"Poisson",sep='.'), paste( colnames(X.g),"Gamma",sep='.'), "\n", sep = "\t")
  }
  fm <- nlminb( start=inits, objective=ldPoisGam.lp, gradient=ldPoisGam.lp.deriv, hessian = NULL, y=y, X.p=X.p, X.g=X.g, offsetty=offset.p, alpha=alpha, control=list(trace=trace), wts=wts1)

  parms <- fm$par
  if( is.null( alpha))
    names( parms) <- c( paste( colnames(X.p),"Poisson",sep='.'), paste( colnames(X.g),"Gamma",sep='.'), "logalpha")
  else 
    names( parms) <- c( paste( colnames(X.p),"Poisson",sep='.'), paste( colnames(X.g),"Gamma",sep='.'))

  if( vcov){
    if( trace!=0)
      print( "Calculating variance matrix of estimates")
    vcovar <- nd2(x0=parms, f=ldPoisGam.lp.deriv, y=y, X.p=X.p, X.g=X.g, offsetty=offset.p, alpha=alpha, wts=wts1)
    vcovar <- 0.5 * ( vcovar + t( vcovar))
    vcovar <- solve( vcovar)
    rownames( vcovar) <- colnames( vcovar) <- names( parms)
  }
  else{
    if( trace!=0)
      print( "Not calculating variance matrix of estimates")
    vcovar <- NULL
  }

  scores <- -ldPoisGam.lp.deriv( parms=parms, y=y, X.p=X.p, X.g=X.g, offsetty=offset.p, alpha=alpha, wts=wts1)

  if( trace !=0)
    print( "Calculating means")
  muLamb <- exp( X.p %*% parms[1:ncol( X.p)] + offset.p)
  muMuZ <- exp( X.g %*% parms[ncol( X.p)+1:ncol( X.g)])
  mu <- muLamb * muMuZ
  fitMu <- cbind( mu, muLamb, muMuZ)
  colnames( fitMu) <- c("total","Lambda","muZ")

  if( residuals){
    if( trace!=0)
      print( "Calculating quantile residuals")
    if( is.null( alpha)){
      resids <- matrix( rep( pPoisGam( y, muLamb, muMuZ, exp( tail( parms,1))), 2), ncol=2)
      resids[y==0,1] <- 0.5 * dPoisGam( y[y==0], muLamb[y==0], muMuZ[y==0], exp( tail( parms,1)), LOG=FALSE)
    }
    else{
      resids <- matrix( rep( pPoisGam( y, muLamb, muMuZ, alpha),2), ncol=2)
      resids[y==0,1] <- 0.5 * dPoisGam( y[y==0], muLamb[y==0], muMuZ[y==0], alpha, LOG=FALSE)
    }
    nzero <- sum( y==0)
    resids[y==0,2] <- runif( nzero, min=rep( 0, nzero), max=2*resids[y==0,1])
    resids <- qnorm( resids)
    colnames( resids) <- c("expect","random")
  }
  else{
    if( trace!=0)
      print( "Not calculating quantile residuals")
    resids <- NULL
  }

  if( trace!=0)
    print( "Done")

  ICadj <- 0
  if( !is.null( alpha))
    ICadj <- ICadj + 1
  AIC <- -2*(-fm$objective) + 2*(length( parms)-ICadj)
  BIC <- -2*(-fm$objective) + log( nrow( X.p))*(length( parms)-ICadj)

  res <- list( coef=parms, logl=-fm$objective, scores=scores, vcov=vcovar, conv=fm$convergence, message=fm$message, niter=fm$iterations, evals=fm$evaluations, call=match.call(), fitted=fitMu, residuals=resids, AIC=AIC, BIC=BIC)

  class( res) <- "pgm"

  return( res)
}


"pPoisGam" <-
function ( q, lambda, mu.Z, alpha) 
{
    tmp <- c( length( q), length( lambda), length(mu.Z), length( alpha))
    names( tmp) <- c( "n","lambda","mu.Z","alpha")
    if( !all( is.element( tmp[-1], c( 1, tmp[1])))) {
      print( "pPoisGam: error -- length of arguments are not compatible")
      return( tmp)
    }

	 #from here on taken from Dunn's function ptweedie.series
	 #I have to admit that I don't quite get it.  I think that it is some sort of quadrature (this is a guess)
	 y <- q
    drop <- 39
    lambdaMax <- N <- max(lambda)
    logfmax <- -log(lambdaMax)/2
    estlogf <- logfmax
    while ((estlogf > (logfmax - drop)) & (N > 1)) {
        N <- max(1, N - 2)
        estlogf <- -lambdaMax + N * (log(lambdaMax) - log(N) + 1) - log(N)/2
    }
    lo.N <- max(1, floor(N))
    lambdaMin <- N <- min(lambda)
    logfmax <- -log(lambdaMin)/2
    estlogf <- logfmax
    while (estlogf > (logfmax - drop)) {
        N <- N + 1
        estlogf <- -lambdaMin + N * (log(lambdaMin) - log(N) + 1) - log(N)/2
    }
    hi.N <- max(ceiling(N))
    cdf <- matrix( 0, nrow=length( y), ncol=1)		#array(dim = length(y), 0)
    for (N in (lo.N:hi.N)) {
        pois.den <- dpois(N, lambda)
        incgamma.den <- pchisq(2 * as.numeric( y) * alpha / as.numeric( mu.Z), 2 * alpha * N)
        cdf <- cdf + pois.den * incgamma.den
    }
    cdf <- cdf + exp(-lambda)
    its <- hi.N - lo.N + 1
 
	return( cdf)

}


"pTweedie" <-
function ( q, mu, phi, p) 
{
  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*mu^(p-1)
  mu.Z <- alpha * tau

  ps <- pPoisGam( q, lambda, mu.Z, alpha)
 
#  require( tweedie)
#  ps <- ptweedie( as.numeric( q), mu=as.numeric( mu), phi=as.numeric( phi), power=as.numeric( p))
  
  return( ps)
}


"rPoisGam" <-
function ( n, lambda, mu.Z, alpha)
{
  mu.N <- lambda
#simulate n random variables from the same compound poisson distribution
  my.fun <- function (parms)
    return( rgamma( n=1, scale=parms[3], shape=parms[1]*parms[2]))
#    return( sum( rgamma( n=parms[1], scale=parms[3], shape=parms[2])))

  tmp <- c( n, length( mu.N), length(mu.Z), length( alpha))
  names( tmp) <- c( "n","mu.N","mu.Z","alpha")
  if( !all( is.element( tmp[-1], c( 1, tmp[1])))) {
    print( "rPoisGam: error -- length of arguments are not compatible")
    return( tmp)
  }
  if( tmp["mu.N"]==1)
    mu.N <- rep( mu.N, tmp["n"])
  if( tmp["mu.Z"]==1)
    mu.Z <- rep( mu.Z, tmp["n"])
  if( tmp["alpha"]==1)
    alpha <- rep( alpha, tmp["n"])

  np <- matrix( rpois( n, mu.N), ncol=1)
  beta <- mu.Z / alpha
  y <- apply( cbind( np, alpha, beta), 1, my.fun)

  return( y)
}


"rTweedie" <-
function ( n, mu, phi, p)
{
  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  alpha <- ( 2-p) / ( p-1)
  tau <- phi*(p-1)*mu^(p-1)
  mu.Z <- alpha * tau
 
  rans <- rPoisGam( n, lambda, mu.Z, alpha)
  return( rans)
}


"simReg" <-
function (n, lambda.tau, mu.Z.tau, alpha, offset1=NULL, X=NULL) 
{
    if (length(lambda.tau) != length(mu.Z.tau) && length(alpha) != 
        1) {
        print("coefficient vectors are inconsistent -- try again")
        return()
    }
    if( is.null( X))
      X <- cbind(1, matrix(rnorm(n * (length(lambda.tau) - 1)), nrow = n))
    else
      X <- as.matrix( X)
    offset2 <- offset1
    if( is.null( offset1))
      offset2 <- rep( 0, nrow( X))
    if( ncol( X) != length( lambda.tau)) {
      print( "Supplied X does not match coefficients")
      X <- cbind(1, matrix(rnorm(n * (length(lambda.tau) - 1)), nrow = n)) }
    lambda <- exp(X %*% lambda.tau + offset2)
    mu <- exp(X %*% mu.Z.tau)
    y <- rPoisGam(n = n, lambda = lambda, mu.Z = mu, alpha = alpha[1])
    if( is.null( offset1)){
      res <- as.data.frame(cbind(y, X))
      colnames(res) <- c("y", "const", paste("x", 1:(length(lambda.tau) - 1), sep = ""))
    }      
    else{
      res <- as.data.frame( cbind( y, offset2, X))
      colnames(res) <- c("y", "offset", "const", paste("x", 1:(length(lambda.tau) - 1), sep = ""))
    }
    attr(res, "coefs") <- list(lambda.tau = lambda.tau, mu.Z.tau = mu.Z.tau, 
        alpha = alpha)
    return(res)
}


"tglm" <-
function ( mean.form, data, wts=NULL, phi=NULL, p=NULL, inits=NULL, vcov=TRUE, residuals=TRUE, trace=1, iter.max=150) 
{
  if( is.null( wts))
    wts <- rep( 1, nrow( data))

  e<-new.env()
  e$wts <- wts
  environment(mean.form) <- e
  temp.p <- model.frame( mean.form, data=as.data.frame( data), weights=wts)
  y <- model.response( temp.p)
  names( y) <- NULL
  X.p <- model.matrix( mean.form, data)
  offset.p <- model.offset( temp.p)
  if ( is.null( offset.p))
    offset.p <- rep( 0, nrow( temp.p))
  wts1 <- as.vector( model.weights( temp.p))

  fm1 <- NULL
  if( is.null( inits)){
    inits <- rep( 0, ncol( X.p))
    if( trace!=0)
      print( "Obtaining initial mean values from log-linear Poisson model -- this might be stupid")
    abit <- 1
    ystar <- round( y, digits=0)
    fm1 <- glm( ystar~-1+X.p, family=poisson( link="log"), weights=wts1)
    inits <- fm1$coef
  }

  if( is.null( phi) & length( inits)==ncol( X.p)){
    if( is.null( fm1))
      inits <- c( inits, 1)
    else{
      if( trace!=0)
        print( "Obtaining initial dispersion from the smaller of the Pearson or Deviance estimator (or 25)")
      if( is.null( p) & length( inits)==ncol( X.p))
        ptemp <- 1.9
      else
        ptemp <- p
      disDev <- fm1$deviance/( length( y) - length( fm1$coef))
      disPear <- sum( ( wts1 * ( y-fm1$fitted)^2) / ( fm1$fitted^ptemp)) / ( length( y) - length( fm1$coef))
      dis <- min( disDev, disPear, 25)
      inits <- c( inits, dis)
    }
  }

  if( is.null( p) & length( inits)==ncol( X.p) + is.null( phi))
    inits <- c( inits, 1.6)

  if( length( inits) != ncol( X.p) + is.null( phi) + is.null( p)) {
    print( "Initial values supplied are of the wrong length -- please check")
    tmp <- c( length( inits), ncol( X.p) + is.null( phi) + is.null( p))
    names( tmp) <- c( "inits", "nParams")
    return( tmp)
  }

  fmTGLM <- tglm.fit( x=X.p, y=y, wts=wts1, offset=offset.p, inits=inits, phi=phi, p=p, vcov=vcov, residuals=residuals, trace=trace, iter.max=iter.max)

  return( fmTGLM)
}


"tglm.fit" <-
function ( x, y, wts=NULL, offset=rep( 0, length( y)), inits=rnorm( ncol( x)), phi=NULL, p=NULL, vcov=TRUE, residuals=TRUE, trace=1, iter.max=150)
{
  if( trace!=0){
    print( "Estimating parameters")
    if( is.null( phi) & is.null( p))
      cat("iter:", "-logl", colnames(x), "phi", "p", "\n", sep = "\t") 
    if( is.null( phi) & !is.null( p))
      cat("iter:", "-logl", colnames(x), "phi", "\n", sep = "\t") 
    if( !is.null( phi) & is.null( p))
      cat("iter:", "-logl", colnames(x), "p", "\n", sep = "\t") 
    if( !is.null( phi) & !is.null( p))
      cat("iter:", "-logl", colnames(x), "\n", sep = "\t") 
  }

  eps <- 1e-5
  my.lower <- rep( -Inf, ncol( x))
  my.upper <- rep( Inf, ncol( x))
  if( is.null( phi))
    {my.lower <- c( my.lower, eps); my.upper <- c( my.upper, Inf)}
  if( is.null( p))
    {my.lower <- c( my.lower, 1+eps); my.upper <- c( my.upper, 2-eps)}
  
  fm <- nlminb( start=inits, objective=ldTweedie.lp, gradient=ldTweedie.lp.deriv, hessian=NULL, lower=my.lower, upper=my.upper, y=y, X.p=x, offsetty=offset, phi=phi, p=p, control=list(trace=trace, iter.max=iter.max), wts=wts)
    
  parms <- fm$par
  tmp <- colnames( x)
  if( !is.null( phi) & !is.null( p))
    tmp <- tmp
  if( !is.null( phi) & is.null( p))
    tmp <- c( tmp, "p")
  if( is.null( phi) & !is.null( p))
    tmp <- c( tmp, "phi")
  if( is.null( phi) & is.null( p))
    tmp <- c( tmp, "phi", "p")
   
  names( parms) <- tmp
  
  if( vcov){
    if( trace!=0)
      print( "Calculating variance matrix of estimates")
    vcovar <- nd2(x0=parms, f=ldTweedie.lp.deriv, y=y, X.p=x, offsetty=offset, phi=phi, p=p, wts=wts)
    vcovar <- 0.5 * ( vcovar + t( vcovar))
    vcovar <- solve( vcovar)
    rownames( vcovar) <- colnames( vcovar) <- names( parms)
  }
  else{
    if( trace!=0)
      print( "Not calculating variance matrix of estimates")
    vcovar <- NULL
  }

  scores <- -ldTweedie.lp.deriv( parms=parms, y=y, X.p=x, offsetty=offset, phi=phi, p=p, wts=wts) 

  if( trace !=0)
    print( "Calculating means")
  mu <- exp( x %*% parms[1:ncol( x)] + offset)

  if( residuals){
    if( trace!=0)
      print( "Calculating quantile residuals")
    if( is.null( phi)){
      phi1 <- parms[ncol( x)+1]
      if( is.null( p))
        p1 <- parms[ncol(x)+2]
      else
        p1 <- p
    }
    else{
      phi1 <- phi
      if( is.null( p))
        p1 <- parms[ncol( x)+1]
      else
        p1 <- p
    }
    resids <- matrix( rep( pTweedie( y, mu, phi1, p1), 2), ncol=2)
    resids[y==0,1] <- 0.5 * dTweedie( y[y==0], mu[y==0], phi1, p1, LOG=FALSE)
    nzero <- sum( y==0)
    resids[y==0,2] <- runif( nzero, min=rep( 0, nzero), max=2*resids[y==0,1])
    resids <- qnorm( resids)
    colnames( resids) <- c("expect","random")
  }
  else{
    if( trace!=0)
      print( "Not calculating quantile residuals")
    resids <- NULL
  }

  if( trace!=0)
    print( "Done")

  ICadj <- 0
  if( !is.null( phi))
    ICadj <- ICadj+1
  if( !is.null( p))
    ICadj <- ICadj+1 

  AIC <- -2*(-fm$objective) + 2*(length( parms)-ICadj)
  BIC <- -2*(-fm$objective) + log( nrow( x))*(length( parms)-ICadj)

  res <- list( coef=parms, logl=-fm$objective, scores=scores, vcov=vcovar, conv=fm$convergence, message=fm$message, niter=fm$iterations, evals=fm$evaluations, call=match.call(), fitted=mu, residuals=resids, AIC=AIC, BIC=BIC)

  class( res) <- "tglm"

  return( res)


}

