#Solves the weighted absolute value norm

dSolver<-function(z, a, extra) {
    x <- z
    if ((is.null(extra$weights)) || (is.null(extra$y))) stop("dSolver needs the additional arguments y and weights!")
    w <- extra$weights
    z <- extra$y
    n <- length(z)
    if (length(a)==0) return(list(x=z,lbd=NULL,f=0,gx=rep(0,n))) 
    if (is.vector(a)) a<-matrix(a,1,2)
    indi<-mkIndi(a,n)
    m<-ncol(indi); h<-rep(0,m)
    for (j in 1:m) {
        ij<-which(indi[,j]==1)
        zj<-z[ij]; wj<-w[ij]
        h[j]<-weightedMedian(zj,wj)
        }
    y<-drop(indi%*%h)
    f<-sum(w*abs(z-y))
    gy<-w*sign(y-z)
    lbd<-mkLagrange(a,gy)
    return(list(x = y, lbd = lbd, f = f, gx = gy))
}