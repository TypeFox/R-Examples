# Author: cns
# basic models

dlmodeler.polynomial <-
function(ord, sigmaH=NA, sigmaQ=0, name=ifelse(ord==0,'level',ifelse(ord==1,'level+trend','polynomial')))
{
	if( ord<0 ) stop("Order must be >= 0")
	m <- ord+1
	if( length(sigmaQ)!=1 & length(sigmaQ)!=m ) stop("SigmaQ has wrong dimension: should be of size ",m)
	d <- 1
	
	a0 <- matrix(0,m,1)
	P0 <- diag(0,m,m)
	P0inf <- diag(m)
	
	Tt <- diag(1,m,m)
	if( m>1 ) for( i in 1:(m-1) ) Tt[i,i+1] <- 1
	Rt <- diag(1,m,m)
	Qt <- diag(sigmaQ^2,m)
	
	Zt <- matrix(c(1,rep(0,m-1)),d,m)
	Ht <- matrix(sigmaH^2,d,d)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



dlmodeler.constant <-
function(cst,name="constant") {
  ret <- dlmodeler.build.polynomial(ord=0,sigmaH=0,sigmaQ=0,name=name)
  ret$a0[1,1] <- cst
  ret$P0inf[1,1] <- 0
  return(ret)
}



dlmodeler.dseasonal <-
		function(ord, sigmaH=NA, sigmaQ=0, name='dseasonal')
{
	if( ord<2 ) stop("Order must be >= 2")
	m <- ord-1
	d <- 1
	
	a0 <- matrix(0,m,1)
	P0 <- diag(0,m,m)
	P0inf <- diag(m)
	
	Tt <- matrix(0,m,m)
	if( m>1 ) for( i in 1:(m-1) ) {
			Tt[1,i] <- -1
			Tt[i+1,i] <- 1
		}
	Tt[1,m] <- -1
	Rt <- matrix(0,m,1)
	Rt[1,1] <- 1
	Qt <- matrix(sigmaQ^2,1,1)
	
	Zt <- matrix(c(1,rep(0,m-1)),d,m)
	Ht <- matrix(sigmaH^2,d,d)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



dlmodeler.tseasonal <-
		function(per, ord=NULL, sigmaH=NA, sigmaQ=0, name='tseasonal')
{
	if( per<=0 ) stop("Period must be > 0")
	if( ((per%%1)!=0) & is.null(ord) ) stop("Order of the trigonometric decomposition must be specified when period is not an integer")
	if( is.null(ord) ) m <- per-1 else m <- 2*ord
	if( (m%%1) != 0 ) stop("Order of the trigonopetric decomposition must be an integer value")
	d <- 1
	
	a0 <- matrix(0,m,1)
	P0 <- diag(0,m,m)
	P0inf <- diag(m)
	
	Tt <- diag(-1,m,m)
	if( per==4 ) {
		M <- matrix( c(0, -1, 1, 0), 2, 2)
	} else {
		f <- 2*base::pi/per
		M <- matrix( c(cos(f),-sin(f),sin(f),cos(f)), 2, 2 )
	}
	N <- M
	if( m>1) for( i in 1:(m/2) ) {
			Tt[(2*i-1):(2*i), (2*i-1):(2*i)] <- N
			N <- M %*% N
		}
	Rt <- diag(1,m,1)
	Qt <- sigmaQ^2*matrix(1,1,1)
	
	Zt <- matrix(rep(c(1,0),m),d,m)
	Ht <- matrix(sigmaH^2,d,d)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



dlmodeler.structural  <-
		function(pol.order=NULL, dseas.order=NULL, tseas.period=NULL, tseas.order=NULL,
				sigmaH=NA, pol.sigmaQ=0, dseas.sigmaQ=0, tseas.sigmaQ=0, name='structural')
{
	if( !is.null(pol.order) ) {
		if( pol.order==0 ) { pol.name <- 'level'
		} else if( pol.order==1 ) { pol.name <- 'level+trend'
		} else pol.name <- 'polynomial'
		mdl1 <- dlmodeler.build.polynomial(pol.order,sigmaH,pol.sigmaQ,name=pol.name)
	} else mdl1 <- NULL
	if( !is.null(dseas.order) ) mdl2 <- dlmodeler.build.dseasonal(dseas.order,0,dseas.sigmaQ,name='seasonal') else mdl2 <- NULL
	if( !is.null(tseas.period) ) mdl3 <- dlmodeler.build.tseasonal(tseas.period,tseas.order,0,tseas.sigmaQ,name='trigonometric') else mdl3 <- NULL
	ret <- mdl1+mdl2+mdl3
  ret$name <- name
  return(ret)
}



dlmodeler.arima <-
		function(ar=c(), ma=c(), d=0, sigmaH=NA, sigmaQ=0, name='arima')
{
	if( d>0 ) stop("case where d>0 is not implemented yet") # TODO
	if( length(ar)==0 & length(ma)==0 ) stop("ar and/or ma terms are missing")
	r <- max(length(ar),length(ma)+1)
	ar <- c(ar,rep(0,r-length(ar)))
	ma <- c(1,ma,rep(0,r-length(ma)-1))
	
	a0 <- matrix(0,r,1)
	P0 <- diag(0,r,r)
	P0inf <- diag(r)
	
	Tt <- cbind(matrix(ar,r,1), diag(1,r,r-1))
	Rt <- matrix(ma,r,1)
	Qt <- matrix(sigmaQ^2,1,1)
	
	Zt <- matrix(c(1,rep(0,r-1)),1,r)
	Ht <- matrix(sigmaH^2,1,1)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



dlmodeler.regression <-
		function(covariates, sigmaH=NA, sigmaQ=0, intercept=FALSE, name='regression')
{
	# covariates must be in hoirontal format (1 row per covariate, as y is formatted)
	if( !is.matrix(covariates) ) covariates <- matrix(covariates,nrow=1)
	n <- NCOL(covariates)
	m <- NROW(covariates)+intercept
	if( intercept ) covariates <- rbind(matrix(1,1,n),covariates)
	if( length(sigmaQ)!=1 & length(sigmaQ)!=m ) stop("SigmaQ has wrong dimension: should be of size ",m)
	d <- 1
	
	a0 <- matrix(0,m,1)
	P0 <- diag(0,m,m)
	P0inf <- diag(m)
	
	Tt <- diag(1,m,m)
	Rt <- diag(1,m,m)
	Qt <- diag(sigmaQ^2,m,m)
	
	Zt <- array(dim=c(d,m,n))
	for( i in 1:n ) Zt[,,i] <- t(covariates[,i])
	Ht <- matrix(sigmaH^2,d,d)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



# classical models
stochastic.level <- function(name="stochastic level") dlmodeler.polynomial(0, sigmaH=0, sigmaQ=NA, name=name)
stochastic.trend <- function(name="stochastic trend") dlmodeler.polynomial(1, sigmaH=0, sigmaQ=c(0,NA), name=name)
stochastic.season <- function(ord, name="stochastic season") dlmodeler.dseasonal(ord, sigmaH=0, sigmaQ=NA, name=name)
random.walk <- function(name="random walk") stochastic.level(name)
deterministic.level <- function(name="deterministic level") dlmodeler.polynomial(0, sigmaH=0, sigmaQ=0, name=name)
deterministic.trend <- function(name="deterministic trend") dlmodeler.polynomial(1, sigmaH=0, sigmaQ=0, name=name)
deterministic.season <- function(ord, name="deterministic season") dlmodeler.dseasonal(ord, sigmaH=0, sigmaQ=0, name=name)



# old function names
dlmodeler.build.polynomial <- dlmodeler.polynomial
dlmodeler.build.constant <- dlmodeler.constant
dlmodeler.build.dseasonal <- dlmodeler.dseasonal
dlmodeler.build.tseasonal <- dlmodeler.tseasonal
dlmodeler.build.structural <- dlmodeler.structural
dlmodeler.build.arima <- dlmodeler.arima
dlmodeler.build.regression <- dlmodeler.regression
