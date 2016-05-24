scio <- function(S, lambda,  thr=1e-4, maxit=1e4, pen.diag=F, sym=T) {
  p <- nrow(S)
  ##w <- matrix(0, p, p)
  w <- diag(diag(S),p)
  nniter <- 0
  jerr <- 0

  isym <- sym*1
  
  lambda <- as.matrix(lambda)
  if (all(nrow(lambda)==p, ncol(lambda)==p)) {
    lambdamat <- lambda
  } else {
    lambdamat <- matrix(lambda[1,1],p,p)
  }

  if (!pen.diag) diag(lambdamat) <- 0.0
  
  mode(p)="integer"
  mode(S)="double"
  mode(w)="double"
  mode(lambdamat)="double"
  mode(thr)="double"
  mode(maxit)="integer"
  mode(isym)="integer"
 
  junk<-.Fortran("scio",
                 p,
                 S,
                 w=w,
                 lambdamat,
                 thr,
                 maxit=maxit,
                 nniter,
                 ierr=integer(1),
                 isym,
                 PACKAGE="scio"
                 )
  
  w = matrix(junk$w, ncol=p)
  return(list(w=w))
}
