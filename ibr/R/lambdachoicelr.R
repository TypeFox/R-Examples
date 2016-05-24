lambdachoicelr <- function(x,ddlobjectif,m=2,s=0,rank,itermax,bs,listvarx) {
  n <- nrow(x)
  d <- ncol(x)
  ddlmin <- choose(m+d-1,m-1)
  if (ddlobjectif<=ddlmin) stop(paste("the objective df is too small, choose it greater than",ddlmin))
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
    vp <- eigen(forwardsolve(t(Rlr),t(forwardsolve(t(Rlr),Slr))),symmetric=TRUE,only.values=TRUE)$values
    ## trace : sum( 1/(1+lambda*vp) )
  trace <- function(loglambda,vp1) {
    sum(1/(1+exp(loglambda)*vp1)) - ddlobjectif
  }
  l1 <- 1
 for (k in 1:25) {
    tr <- sum( 1/(1+l1*vp) )
       if (tr <= ddlobjectif)  break
        l1 <- l1 * 4
    }
     l2 <- 1
    for (k in 1:25) {
        tr <- sum( 1/(1+l2*vp) )
         if (tr > ddlobjectif) break
        l2 <- l2/4
    }
    resultat <- uniroot(trace,c(log(l2),log(l1)),vp1=vp,maxiter =itermax)
    return(exp(resultat$root))
}
