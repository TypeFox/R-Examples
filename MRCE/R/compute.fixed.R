compute.fixed=function(X,Y, lam2, omega, tol.in, maxit.in, silent)
{
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  if(!is.matrix(lam2))
    nlam=matrix(n*lam2, nrow=p, ncol=q) else nlam=n*lam2  
  mx=apply(X, 2, mean)
  my=apply(Y, 2, mean)
  X=scale(X, center=mx, scale=FALSE)
  Y=scale(Y, center=my, scale=FALSE)
  xty=crossprod(X,Y)
  xtx=crossprod(X)
  tolmult=sum(crossprod(Y)*omega)
  tol.in=tol.in*tolmult
  xtyom=xty%*%omega
  new.B=rblasso(s=xtx, m=xtyom, om=omega, nlam=nlam, n=n, B0=NULL, soft=NULL, tol=tol.in, maxit=maxit.in, quiet=silent)$B
  muhat=as.numeric(my - crossprod(new.B, mx))
  return(list(Bhat=new.B, muhat=muhat, omega=omega, mx=mx, my=my))
}

