#Chebyshev norm

mSolver<-function(z, a, extra) {
    x <- z
    if ((is.null(extra$weights)) || (is.null(extra$y))) stop("mSolver needs the additional arguments y and weights!")
    w<-extra$weights
    z<-extra$y
    n<-length(z)
    if (length(a)==0) return(list(x=z,lbd=NULL,f=0,gx=rep(0,n))) 
    if (is.vector(a)) a<-matrix(a,1,2)
    indi<-mkIndi(a,n)
    m<-ncol(indi); h<-rep(0,m)
    for (j in 1:m) {
        ij<-which(indi[,j]==1)
        zj<-z[ij]; wj<-w[ij]
        h[j]<-weightedMidRange(zj,wj)
        }
    y<-drop(indi%*%h); dv<-w*(y-z)
    i1<-which.max(dv); i2<-which.min(dv)
    f<-max(abs(dv))
    gy1<-rep(0,n); gy1[i1]<-w[i1]
    lbd1<-mkLagrange(a,gy1)
    gy2<-rep(0,n); gy2[i2]<--w[i2]
    lbd2<-mkLagrange(a,gy2)
    lbd<-(w[i2]*lbd1+w[i1]*lbd2)/(w[i1]+w[i2])
    gy<-(w[i2]*gy1+w[i1]*gy2)/(w[i1]+w[i2])
    return(list(x = y, lbd = lbd, f = f, gx = gy))
}