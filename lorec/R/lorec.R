lorec <- function(Sig, L=NULL, S=NULL, lambda, delta,  thr=1.0e-4, maxit=1e4) {

  p <- nrow(Sig)
  mode(p)="integer"
  if (is.null(L)) L <- diag(diag(Sig))
  if (is.null(S)) S <- diag(diag(Sig))
  
  mode(L)="double"
  mode(S)="double"
  
  mode(Sig)="double"  
  mode(lambda)="double"
  mode(delta)="double"
  mode(thr)="double"
  mode(maxit)="integer"

  junk<-.Fortran("lorec",
                 p,
                 Sig,
                 L=L,
                 S=S,
                 lambda,
                 delta,
                 thr,
                 maxit=maxit,
                 ierr=integer(1),
                 PACKAGE="lorec"
                 )
  
  L = matrix(junk$L, ncol=p)
  S= matrix(junk$S, ncol=p)
  return(list(L=L,S=S))
}
