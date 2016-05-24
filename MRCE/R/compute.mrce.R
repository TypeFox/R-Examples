compute.mrce=function(X,Y, lam1, lam2, tol.out, tol.in, maxit.out, maxit.in, silent,
                      cov.tol, cov.maxit,informed=NULL, eps)
{
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  if(is.null(informed))
  {
    if(!is.matrix(lam2))
      nlam=matrix(n*lam2, nrow=p, ncol=q) else nlam=n*lam2
    
    mx=apply(X, 2, mean)
    my=apply(Y, 2, mean)
    X=scale(X, center=mx, scale=FALSE)
    Y=scale(Y, center=my, scale=FALSE)
    yty=crossprod(Y)
    xty=crossprod(X,Y)
    xtx=crossprod(X)
    old.B=matrix(0, nrow=p, ncol=q)
    tolmult=sum(diag(yty)/n)
    tout=tol.out*tolmult
    residual.cov = yty/n
    sigma=diag(diag(residual.cov))
    om=diag(1/diag(residual.cov))
    omoff=om
    diag(omoff)=0
    old.obj=sum(residual.cov*om)-determinant(om, logarithm=TRUE)$mod[1] + lam1*sum(abs(omoff))+2*sum(lam2*abs(old.B))
  } else
  {
    nlam=matrix(n*lam2, nrow=p, ncol=q)
    mx=informed$mx
    my=informed$my
    xtx=informed$xtx
    xty=informed$xty
    yty=informed$yty
    old.B=informed$Bhat
    om=informed$omega
    sigma=informed$sigma
    tolmult=sum(diag(yty)/n)
    tout=tol.out*tolmult
    residual.cov = crossprod(Y-X%*%old.B)/n
    omoff=om
    diag(omoff)=0
    old.obj=sum(residual.cov*om)-determinant(om, logarithm=TRUE)$mod[1] + lam1*sum(abs(omoff))+2*sum(lam2*abs(old.B))
  }
  Moff=lam1*(1-diag(q))
  k=0  
  iterating=TRUE
  while(iterating)
  {
    k=k+1  
    if(min(diag(residual.cov)) < eps)
    {
      cat("A perfect fit occured (lam2 may be too small). Terminated early.\n")
      break
    }   
    if(k == 1)
    {
      cov.out=NULL
      cov.out$X=om
      cov.out$W=sigma
    }
    cov.out=QUIC(S=residual.cov, rho=Moff,tol=cov.tol, msg=0, maxIter=cov.maxit, X.init=cov.out$X, W.init=cov.out$W)
    if(!is.matrix(cov.out$X))
      cov.out$X=matrix(cov.out$X, nrow=q, ncol=q)
    if(!is.matrix(cov.out$W))
      cov.out$W=matrix(cov.out$W, nrow=q, ncol=q)  
    om=cov.out$X   
    tolinmult=sum(yty*om)/n 
    ## added
    omoff=om
    diag(omoff)=0
    obj.after.omega=sum(residual.cov*om)-determinant(om, logarithm=TRUE)$mod[1] + lam1*sum(abs(omoff))+2*sum(lam2*abs(old.B))
    if(!silent) cat("k =", k, "obj. fn. val. after Omega update is", obj.after.omega, "\n")  
	xtyom=xty%*%om  
    soft=xtyom - xtx%*%old.B%*%om + old.B*tcrossprod(diag(xtx), diag(om)) 
    outlasso=rblasso(s=xtx, m=xtyom, om=om, nlam=nlam, n=n,B0=old.B, soft=soft, objective=obj.after.omega, tol=(tol.in*tolinmult), maxit=maxit.in, quiet=silent)		
    old.B=outlasso$B
    residual.cov = crossprod(Y-X%*%old.B)/n
    new.obj=outlasso$obj
    bdist = old.obj-new.obj
    iterating= (bdist > tout) & (k <= maxit.out)
    old.obj=new.obj
    if(!silent) cat("k =", k, "obj. fn. val. after B update is", new.obj, "\n")
  }
  if(!silent) cat("Total outer iterations for MRCE : ", k, "\n")
  muhat=as.numeric(my - crossprod(old.B, mx))
  return(list(Bhat=old.B, muhat=muhat, omega=om, sigma=cov.out$W, mx=mx, my=my))
}

