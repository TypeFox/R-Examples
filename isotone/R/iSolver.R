# SILF Loss

# SILF Loss (Chu et al., 2004)

iSolver<-function(z, a, extra) {
  x <- z
    if ((is.null(extra$weights)) || (is.null(extra$y)) || (is.null(extra$eps)) || (is.null(extra$beta))) stop("iSolver needs the additional arguments y, weights, eps, and beta!")
 
  w <- extra$weights
  z <- extra$y
  eps <- extra$eps
  beta<-extra$beta
  fobj<-function(x) {
    delta<-abs(x-z)
    g<-((delta-(1-beta)*eps)^2)/(4*beta*eps)
    g[which(delta < (1-beta)*eps)]<-0
    ii<-which(delta > (1+beta)*eps)
    g[ii]<-delta[ii]-eps
    return(sum(w*g))
    }
  gobj<-function(x) {
    y<-x-z
    g<-rep(0,length(y))
    g[which(y < -(1+beta)*eps)]<--1
    ii<-which((y > -(1+beta)*eps) & (y < -(1-beta)*eps))
    g[ii]<-(y[ii]+(1-beta)*eps)/(2*beta*eps)
    ii<-which((y > (1-beta)*eps) & (y < (1+beta)*eps))
    g[ii]<-(y[ii]-(1-beta)*eps)/(2*beta*eps)
    g[which(y > (1+beta)*eps)]<-1
    return(w*g)
    }
return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}