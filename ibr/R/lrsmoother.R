lrsmoother <- function(x,bs,listvarx,lambda,m,s,rank) {
    if (bs=="tp") { 
        objet <- list(term=listvarx,bs.dim=rank,fixed=FALSE,dim=ncol(x),p.order=m,by="NA",label=paste("s(",paste(listvarx,collapse=","),")",sep=""),xt=NULL,id=NULL,sp=NULL)
        attr(objet,"class") <- "tp.smooth.spec"
    } else {
        objet <- list(term=colnames(x),bs.dim=rank,fixed=FALSE,dim=ncol(x),p.order=c(m,s),by="NA",label=paste("s(",paste(listvarx,collapse=","),")",sep=""),xt=NULL,id=NULL,sp=NULL)
        attr(objet,"class") <- "ds.smooth.spec"       
    }
    sm <- smoothCon(objet,data=data.frame(x),knots=NULL)[[1]]
    Xlr <- sm$X
    Slr <- sm$S[[1]]
    qrx <- qr(Xlr)
    Rlr <- qr.R(qrx)
    es <- eigen(forwardsolve(t(Rlr),t(forwardsolve(t(Rlr),Slr))),symmetric=TRUE)
    ## dans l'ordre decroissant 
    eigenvaluesS1 <- rev(1/(1+lambda*es$values))
    U <- qr.Q(qrx)%*%es$vectors[,rank:1]
    ## calcul de (R^-1)U
    Rm1U <- backsolve(qr.R(qrx),es$vectors[,rank:1])
    return(list(vectors=U,values=eigenvaluesS1,Rm1U=Rm1U,smoothobject=sm))
}
