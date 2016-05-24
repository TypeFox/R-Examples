tune.cv=function(X, F, d, lam.vec, kfold=5, silent=TRUE, qrtol=1e-16, cov.tol=1e-4, cov.maxit=1e3)
{  
  ## X ~ F, both uncentered
  n=nrow(X)
  p=ncol(X)
  r=ncol(F)
  ind = sample(n)
  err.mat = array(0, length(lam.vec))
  for (k in 1:kfold)
  {
    foldind = ind[ (1+floor((k-1)*n/kfold)):floor(k*n/kfold) ]
    F.tr=F[-foldind, , drop=FALSE]
    X.tr=X[-foldind, , drop=FALSE]
    F.va=F[foldind, , drop=FALSE]
    X.va=X[foldind, , drop=FALSE]
    mtrf=apply(F.tr, 2, mean)
    F.tr=scale(F.tr, scale=FALSE, center=mtrf)
    F.va=scale(F.va, scale=FALSE, center=mtrf)
    mtrx=apply(X.tr, 2, mean)
    X.tr=scale(X.tr, scale=FALSE, center=mtrx)
    X.va=scale(X.va, scale=FALSE, center=mtrx)
    n.tr=dim(X.tr)[1]
    
    FtF.tr=crossprod(F.tr)
    Bhat.tr = qr.solve(FtF.tr, crossprod(F.tr, X.tr), tol=qrtol)
    S.tr=crossprod(X.tr - F.tr%*%Bhat.tr)/n.tr
    sdi.tr=1/(sqrt(diag(S.tr)))
    R.tr=sdi.tr * S.tr * rep(sdi.tr, each = p)
    Phi.tr=FtF.tr/n.tr
    if( r > 1 )
    {
      e.out=eigen(Phi.tr, symmetric=TRUE)
      Phi.tr.sq =tcrossprod(e.out$vec*rep(e.out$val^(0.5), each=r),e.out$vec)
      Phi.tr.nsq =tcrossprod(e.out$vec*rep(e.out$val^(-0.5), each=r),e.out$vec)
    }else
    {
      Phi.tr.sq = Phi.tr^(0.5)
      Phi.tr.nsq= Phi.tr^(-0.5)
    }
    U=Phi.tr.sq%*%Bhat.tr
    for(i in 1:length(lam.vec))
    {
      lam=lam.vec[i]
      Moff=matrix(lam, nrow=p, ncol=p)
      diag(Moff)=0      
      if(i == 1)
      {
        cov.out=NULL
        cov.out$X=diag(p)
        cov.out$W=diag(p)
      }
      if(lam < 1 )
      {  
        cov.out=QUIC(S=R.tr, rho=Moff, tol=cov.tol, msg=(1*(!silent)), maxIter=cov.maxit, X.init=cov.out$X, W.init=cov.out$W)
      }  else
      {
        cov.out$X=diag(p)
        cov.out$W=diag(p)
      }
      if(!is.matrix(cov.out$X))
        cov.out$X=matrix(cov.out$X, nrow=p, ncol=p)
      if(!is.matrix(cov.out$W))
        cov.out$W=matrix(cov.out$W, nrow=p, ncol=p)  
      W=sdi.tr * cov.out$X * rep(sdi.tr, each = p)      

      K = tcrossprod(U%*%W, U)
      Vd = eigen(K, symmetric=TRUE)$vec[, 1:d, drop=FALSE]
      B = crossprod(Vd, Phi.tr.nsq)
      G = crossprod(U, Vd)
      val.cov=crossprod(X.va-tcrossprod(F.va, G%*%B))/(dim(X.va)[1]) 
      out.det=determinant(W, logarithm=TRUE)
      logdet=out.det$mod[1]
      if(out.det$sign > 0)
      {
        err.mat[i]=err.mat[i] + sum(val.cov*W) - logdet
      }else
      {
        err.mat[i]=err.mat[i]+Inf
      }			
    }
    if(!silent) cat("finished fold k = ", k, "\n")
  }
  best.lam = lam.vec[which.min(err.mat)]
  return(list(best.lam=best.lam, err.vec=err.mat))
} 
