#solver the general LS norm (y-x)'W(y-x)

lfSolver<-function(z, a, extra) {
    x <- z
    if ((is.null(extra$weights)) || (is.null(extra$y))) stop("lfSolver needs the additional arguments y and weights!")
    w <- extra$weights
    z <- extra$y
    n <- length(z)
    if (length(a)==0) return(list(x=z,lbd=NULL,f=0,gx=rep(0,n))) 
    if (is.vector(a)) a<-matrix(a,1,length(a))
    indi<-mkIndi(a,n)
    h<-crossprod(indi,w%*%indi); r<-drop(crossprod(indi,w%*%z))
    b<-solve(h,r); y<-drop(indi%*%b); gy<-2*drop(w%*%(y-z))
    lbd<-mkLagrange(a,gy)
    f<-sum(w*outer(y-z,y-z))
    return(list(x = y, lbd = lbd, f = f, gx = gy))
}