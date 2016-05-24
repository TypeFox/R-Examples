mrce=function(X,Y, lam1=NULL, lam2=NULL, lam1.vec=NULL, lam2.vec=NULL,
     method=c("single", "cv", "fixed.omega"),
     cov.tol=1e-4, cov.maxit=1e3, omega=NULL, 
     maxit.out=1e3, maxit.in=1e3, tol.out=1e-8, 
     tol.in=1e-8, kfold=5, silent=TRUE, eps=1e-5)
{
  method=match.arg(method)
  out=switch(method, single=compute.mrce(X=X,Y=Y, lam1=lam1, lam2=lam2, tol.out=tol.out, 
                                        tol.in=tol.in, maxit.out=maxit.out, maxit.in=maxit.in, 
                                        silent=silent, cov.tol=cov.tol, 
                                        cov.maxit=cov.maxit, informed=NULL, eps=eps),
                     fixed.omega=compute.fixed(X=X,Y=Y, lam2=lam2, omega=omega, tol.in=tol.in,
                                        maxit.in=maxit.in, silent=silent),
                     cv=mrce.cv(X=X, Y=Y, lam.vec.1=lam1.vec, lam.vec.2=lam2.vec, 
                                      kfold=kfold,tol.out=tol.out, tol.in=tol.in, maxit.out=maxit.out, 
                                      maxit.in=maxit.in, silent=silent,
                                      cov.tol=cov.tol, cov.maxit=cov.maxit, eps=eps))
  return(out)  
}

