Ssimpute.svd <-
  function (x, J = 2, thresh = 1e-05,lambda=0,maxit=100,trace.it=FALSE,warm.start=NULL,...) 
{
###This function expects an object of class "Incomplete" which inherits from "sparseMatrix", where the missing entries
###are replaced with zeros. If it was centered, then it carries the centering info with it
  this.call=match.call()
  a=names(attributes(x))
  binames=c("biScale:row","biScale:column")
  if(all(match(binames,a,FALSE))){
    biats=attributes(x)[binames]
  } else biats=NULL
  
  if(!inherits(x,"dgCMatrix"))x=as(x,"dgCMatrix")
  irow=x@i
  pcol=x@p
  n <- dim(x)
  m <- n[2]
  n <- n[1]
  nz=nnzero(x)
  xres=x
   if(!is.null(warm.start)){
###must have u,d and v components
    if(!all(match(c("u","d","v"),names(warm.start),0)>0))stop("warm.start does not have components u, d and v")
    D=warm.start$d
    nzD=sum(D>0)
    JD=min(nzD,J)
    U=warm.start$u[,seq(JD),drop=FALSE]
    V=warm.start$v[,seq(JD),drop=FALSE]
    D=D[seq(JD)]
    BD=UD(V,D,m)
    xhat.Omega=suvC(U,BD,irow,pcol)
    xres@x=x@x-xhat.Omega
    xfill=splr(xres,U,BD)
    svd.xfill=Ssvd.als(xfill,J,thresh,lambda,maxit,trace.it,warm.start=list(u=U,d=D,v=V))
    }
  else  {
    xfill=x
    svd.xfill=Ssvd.als(xfill,J,thresh,lambda,maxit,trace.it)
  }
  ratio <- 1
  iter <- 0
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    svd.old=svd.xfill
    BD=UD(svd.xfill$v,svd.xfill$d,m)
    U=svd.xfill$u
    xhat.Omega=suvC(U,BD,irow,pcol)
    xres@x=x@x-xhat.Omega
    xfill=splr(xres,U,BD)
    svd.xfill=Ssvd.als(xfill,J,thresh,lambda,maxit,trace.it,warm.start=svd.xfill)
    ratio=Frob(svd.old$u,svd.old$d,svd.old$v,svd.xfill$u,svd.xfill$d,svd.xfill$v)
    if(trace.it){
      obj=(.5*sum(xres@x^2)+lambda*sum(svd.old$d))/nz
    cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
    }
  }
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
   ###Final cleanup of svd
    A=xfill%*%svd.xfill$v
    Asvd=svd(A)
    U=Asvd$u
    Dsq=Asvd$d
    V=svd.xfill$v %*% Asvd$v
    Dsq=pmax(Dsq-lambda,0)
 ###
   J=min(sum(Dsq>0)+1,J)
  svd.xfill=list(u=U[, seq(J)], d=Dsq[seq(J)], v=V[,seq(J)])
    attributes(svd.xfill)=c(attributes(svd.xfill),list(lambda=lambda,call=this.call),biats)
  svd.xfill
}
