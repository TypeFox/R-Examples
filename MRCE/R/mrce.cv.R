mrce.cv=function(X, Y, lam.vec.1, lam.vec.2, kfold=5,
        tol.out, tol.in, maxit.out, maxit.in, silent,
        cov.tol, cov.maxit, eps)
{
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  err.mat = matrix(0, nrow=length(lam.vec.1), ncol=length(lam.vec.2))
  ind=sample(n)
  for (k in 1:kfold)
  {
    if(!silent) cat("Staring fold : ", k, "\n")
    foldind = ind[ (1+floor((k-1)*n/kfold)):floor(k*n/kfold) ]
    X.tr=X[-foldind,]
    Y.tr=Y[-foldind,]
	X.va=X[foldind,,drop=FALSE]
	Y.va=Y[foldind,,drop=FALSE ]
 
    mtrx=apply(X.tr, 2, mean)
    X.tr=scale(X.tr, scale=FALSE, center=mtrx)
    X.va=scale(X.va, scale=FALSE, center=mtrx)
    mtry=apply(Y.tr, 2, mean)
    Y.tr=scale(Y.tr, scale=FALSE, center=mtry)
    Y.va=scale(Y.va, scale=FALSE, center=mtry)
    
    n.tr=nrow(Y.tr)
    informed=NULL
    informed$mx=mtrx
    informed$my=mtry
    informed$yty=crossprod(Y.tr)
    informed$xtx=crossprod(X.tr)
    informed$xty=crossprod(X.tr, Y.tr)
    informed$Bhat=matrix(0, nrow=p, ncol=q)
    informed$omega=diag(n.tr/diag(informed$yty))
    informed$sigma=diag(diag(informed$yty)/n.tr)
    
	for(i in 1:length(lam.vec.1)) 
    {
      if(i > 1)
        informed=first
      for(j in 1:length(lam.vec.2))
      {
      tmp.out=compute.mrce(X=X.tr,Y=Y.tr, lam1=lam.vec.1[i], lam2=lam.vec.2[j], 
                           tol.out=tol.out, tol.in=tol.in, maxit.out=maxit.out, 
                           maxit.in=maxit.in, silent=silent,
                           cov.tol=cov.tol, cov.maxit=cov.maxit, 
                           informed=informed, eps=eps)
        err.mat[i,j]=err.mat[i,j]+mean((Y.va-X.va%*%tmp.out$Bhat)^2)
        informed$Bhat=tmp.out$Bhat
        informed$omega=tmp.out$omega
        informed$sigma=tmp.out$sigma  
        if(j==1)
          first=informed        
      }
    }
  }
  ## find the (i,j) for the minimum of err.mat
  tmp = which.min(err.mat) %% (dim(err.mat)[1])
  tmp = (tmp != 0)*tmp + (tmp == 0)*(dim(err.mat)[1])
  best.i=tmp
  best.j=which.min(err.mat[tmp,])
  best.lam1=lam.vec.1[best.i]
  best.lam2=lam.vec.2[best.j]
  
  out=compute.mrce(X=X,Y=Y, lam1=best.lam1, lam2=best.lam2, 
                           tol.out=tol.out, tol.in=tol.in, maxit.out=maxit.out, 
                           maxit.in=maxit.in, silent=silent,
                           cov.tol=cov.tol, cov.maxit=cov.maxit, 
                           informed=NULL, eps=eps)
  return(list(Bhat=out$Bhat, muhat=out$muhat, omega=out$omega, mx=out$mx, my=out$my, 
         best.lam1=best.lam1, best.lam2=best.lam2, cv.err=err.mat))
}

