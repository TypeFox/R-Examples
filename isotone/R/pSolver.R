#Quantile Regression

pSolver<-function(z, a, extra) {
    x <- z
     if ((is.null(extra$weights)) || (is.null(extra$y)) || (is.null(extra$aw)) || (is.null(extra$bw))) stop("pSolver needs the additional arguments y, weights, aw, bw!")
    w<-extra$weights
    z<-extra$y
    aw <- extra$aw
    bw <- extra$bw
    n<-length(z)
    if (length(a)==0) return(list(x=z,lbd=NULL,f=0,gx=rep(0,n))) 
    if (is.vector(a)) a<-matrix(a,1,2)
    indi<-mkIndi(a,n)
    m<-ncol(indi); h<-rep(0,m)
    for (j in 1:m) {
        ij<-which(indi[,j]==1)
        zj<-z[ij]; wj<-w[ij]
        h[j]<-weightedFractile(zj,wj,aw,bw)
        }
    y<-drop(indi%*%h); dv<-ifelse(y<=z,w*aw*(z-y),w*bw*(y-z))
    f<-sum(dv); gy<-ifelse(y<=z,-w*aw,w*bw)
    lbd<-mkLagrange(a,gy)
    return(list(x = y, lbd = lbd, f = f, gx = gy))
}