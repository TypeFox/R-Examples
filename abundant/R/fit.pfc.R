fit.pfc=function(X, y, r=4, d=NULL, F.user=NULL, weight.type=c("sample", "diag", "L1"), 
        lam.vec=NULL, kfold=5, silent=TRUE, qrtol=1e-10, cov.tol=1e-4, cov.maxit=1e3, NPERM=1e3, level=0.01)
{
  weight.matrix=match.arg(weight.type)
  X=as.matrix(X)
  y=as.numeric(y)
  p=ncol(X)
  n=nrow(X)
  WML1=weight.matrix=="L1"
  if( is.null(d) & WML1 )
    stop("d must be specified when the L1 weight matrix is used\n")
  if(is.null(lam.vec) & WML1)
    stop("lam.vec must be specified when the L1 weight matrix is used\n")

  if(is.null(F.user))  ## use polynomial basis functions
  {
    f=NULL;
    for(k in 1:r) f=cbind(f, y^k);
  } else
  {
    f=F.user
    r=ncol(f);  
  }
  if(!is.null(d))
    if(d > r)
      stop("d cannot exceed r\n")

  if(n < (r+2))
    stop("that sample size n must be greater than or equal to r+2\n")

  best.lam=NULL
  err.vec=NULL
  if(!is.null(lam.vec))
    lam.vec=sort(lam.vec,  decreasing=TRUE)
  if(weight.matrix=="L1")
  {  
    if(length(lam.vec) > 1 )
    {
      cvout=tune.cv(X=X,F=f, d=d, lam.vec=lam.vec, kfold=kfold, silent=silent, qrtol=qrtol, cov.tol=cov.tol, cov.maxit=cov.maxit)
      best.lam=cvout$best.lam
      err.vec=cvout$err.vec
    } else
    {
      best.lam=lam.vec
    }
  } 
  mx=apply(X, 2, mean)
  mf=apply(f, 2, mean)
  Xc=scale(X, center=mx, scale=FALSE)
  fc=scale(f, center=mf, scale=FALSE)
  ftf=crossprod(fc)
  Bhat = qr.solve(ftf, crossprod(fc, Xc), tol=qrtol)
  res=Xc-fc%*%Bhat
  nt=n-r-1
  Deltahat=crossprod(res)/nt
  if(weight.matrix=="sample")
  {
    if(p > nt) 
    {    
      svdout=try(expr=svd(res), silent=TRUE)
      if( inherits(svdout, "try-error") )
      { 
        ## using the eigendecomposition
        eiout=eigen(Deltahat, symmetric=TRUE)
        evs=c(eiout$val[1:nt]^(-1), rep(0, p-nt))
        W=tcrossprod(eiout$vec * rep(evs, each=p), eiout$vec)
      } else
      {
        evs=c(svdout$d[1:nt]^(-2), rep(0, min(c(n,p))-nt))
        W = nt*tcrossprod(svdout$v*rep(evs, each=p),svdout$v)
      }
    } else
    {
      W=qr.solve(Deltahat, tol=qrtol)
    }
  } else if (weight.matrix=="diag")
  {
    W=diag(1/diag(Deltahat))
  } else  ## L1
  {
    sdi=1/(sqrt(diag(Deltahat)))
    Corhat=sdi * Deltahat * rep(sdi, each = p)
    Moff=matrix(best.lam, nrow=p, ncol=p)
    diag(Moff)=0
    cov.out=QUIC(S=Corhat, rho=Moff, tol=cov.tol, msg=(1*(!silent)), maxIter=cov.maxit)
    W=cov.out$X
    if(!is.matrix(W)) W=matrix(W, nrow=p, ncol=p)
    W=sdi * W * rep(sdi, each = p)   
  }
  Phi=ftf/n
  test.info=NULL
  if(is.null(d))
  {
    ## select d
    stat.list=NULL
    d0.list=NULL
    pv.list=NULL
    for(k in 1:r)
    {
      d0=k-1
      d=k
      test=ptest(X=X, Xc=Xc,y=y, f=f, Bhat=Bhat,weight.matrix=weight.matrix, fc=fc, Phi=Phi, W=W, d0=d0, NPERM=NPERM, qrtol=qrtol)
      pv=mean(test$ndist <= test$stat)
      d0.list[k]=d0
      pv.list[k]=pv
      stat.list[k]=test$stat
      if(pv > level)
      {
        d=d0      
        break
      }        
    }
    test.info=data.frame(d0=d0.list, test.statistic=stat.list, pvalue=pv.list)   
  }
  ## compute the reduced rank solution
  if(d > 0)
  {
    out=solve.exact(Phi=Phi, Bhat=Bhat, d=d, W=W, qrtol=qrtol)
  } else
  {
    out=NULL
  }
  return(list(Gamhat=out$G, bhat=out$B, Rmat=out$Rmat, What=W, d=d, r=r, GWG=out$GWG, fc=fc, Xc=Xc, y=y, 
              mx=mx, mf=mf, best.lam=best.lam, lam.vec=lam.vec, err.vec=err.vec, test.info=test.info))
} 
