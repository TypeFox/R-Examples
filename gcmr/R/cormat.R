## cormat must have class cormat.gcmr and the following elements:
## - npar: number of parameters.
## - start(): initial estimates; optional attributes lower and upper can be used
##   to specify box constrained parameters
## - chol(tau,not.na): returns the cholesky factor of the correlation matrix
##   (only for the not.na observations; not.na is a logical vector)
##   This function should return NULL if tau is outside parameter space

## working independence correlation
ind.cormat <- function() {
    ans <- list()
    ans$npar <- 0
    ans$independent <- TRUE
    ans$start <- function() double(0)
    ans$chol <- function(tau , not.na) diag(rep(1,sum(not.na)))
    class( ans ) <- c("ind.gcmr", "cormat.gcmr")
    ans
}

## arma(p,q) correlation for time-series
arma.cormat <- function( p=0 , q=0 ) {
    if(p==0 && q==0)
      return( ind.cormat() )
    iar <- if ( p ) 1:p else NULL
    ima <- if ( q ) (p+1):(p+q) else NULL
    ans <- list()
    ans$npar <- p+q
    ans$start <- function() {
        tau <- rep(0, p+q)
        names(tau) <- c(if ( p ) paste("ar",1:p,sep="") else NULL ,
                        if ( q ) paste("ma",1:q,sep="") else NULL )
        tau
    }
    ans$chol <- function( tau , not.na ) {
        if ( ( p && any(Mod(polyroot(c(1,-tau[iar])))<1.01) ) ||
            ( q && any(Mod(polyroot(c(1, tau[ima])))<1.01) ) )
            return( NULL )
        n <- length(not.na)
        rho <- ARMAacf(tau[iar],tau[ima],n-1)
        r <- seq(1,n)[not.na]
        chol(outer( r , r , function(i,j) rho[1+abs(i-j)] ))
    }
    class( ans ) <- c("arma.gcmr", "cormat.gcmr")
    ans
}

## clustered data
## assume that it is not possible that all the observations inside a cluster
## can be missing
cluster.cormat <- function(id, type=c("independence", "ar1", "ma1",
                                      "exchangeable", "unstructured")) {
    type <- match.arg(type)
    if(!length(rle(id)$values)==length(unique(id)))
        stop("data must be sorted in way that observations from the same cluster are contiguous")
    ng <- 1:length(unique(id))
    if (!(length(ng)>1)) stop("only one strata")
    if (type=="independence") {
        ans <- ind.cormat()
        ans$id <- id
        return(ans)
    }
    ans <- list(type=type,id=id)
    ans$npar <-  if(type!="unstructured") 1 else choose(max(table(id)), 2) 
    data <- data.frame(id=id)
    fn <- switch(type,
                 "ar1"=function(g) nlme::corAR1(g, form= ~1|id),
                 "ma1"=function(g) nlme::corARMA(g, form= ~1|id, p=0, q=1),
                 "exchangeable"=function(g) nlme::corCompSymm(g, form= ~1|id),
                 "unstructured"=function(g) nlme::corSymm(g, form= ~1|id))
    ans$start <- function() {
        np <-  if(type!="unstructured") 1 else choose(max(table(id)), 2) 
        tau <- rep(0, np)
        names(tau) <- switch(type, "ar1"="ar1", "ma1"="ma1", "exchangeable"="tau",
                             "unstructured"=paste("tau", 1:ans$npar, sep=""))
        eps <- sqrt(.Machine$double.eps)
        attr(tau,"lower") <- rep(-1+eps,np)
        attr(tau,"upper") <- rep(1-eps,np)
        tau
    }
    ans$chol <- function(tau, not.na) {
        q <- try(nlme::corMatrix(nlme::Initialize(fn(tau),data=data)),silent=TRUE)
        if (inherits(q,"try-error")) return(NULL)
        g <- split(not.na,id)
        q <- try(lapply(ng,function(i) chol(q[[i]][g[[i]],g[[i]]])),silent=TRUE)
        if (inherits(q,"try-error") ) NULL else q
    }
    class( ans ) <- c("cluster.gcmr", "cormat.gcmr")
    ans
}

## Matern correlation for spatial data
## D is a distance matrix, alpha is the smoothing parameter
matern.cormat <- function(D, alpha=0.5) {
    ans <- list()
    ans$npar <- 1
    ans$start <- function() {
        tau <- median(D)
        names(tau) <- c("tau")
        attr(tau,"lower") <- sqrt(.Machine$double.eps)
        tau
    }
    ans$chol <- function( tau, not.na ){
        S <- geoR::matern(D, tau, alpha)
        q <- try(chol(S[not.na,not.na]),silent=TRUE)
        if( inherits(q,"try-error") ) NULL else q
    }
    class( ans ) <- c("matern.gcmr", "cormat.gcmr")
    ans
}
