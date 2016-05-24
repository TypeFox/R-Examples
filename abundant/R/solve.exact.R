solve.exact = function(Phi, Bhat, d, W, qrtol)
{
  r=nrow(Phi)
  if( r  > 1 )
  {
    e.out=eigen(Phi, symmetric=TRUE)
    Phi.sq =tcrossprod(e.out$vec*rep(e.out$val^(0.5), each=r),e.out$vec)
    Phi.nsq =tcrossprod(e.out$vec*rep(e.out$val^(-0.5), each=r),e.out$vec)
  }else
  {
    Phi.sq = Phi^(0.5)
    Phi.nsq= Phi^(-0.5)
  }
  U=Phi.sq%*%Bhat
  K = tcrossprod(U%*%W, U)
  Vd = eigen(K, symmetric=TRUE)$vec[, 1:d, drop=FALSE]
  B = crossprod(Vd, Phi.nsq)
  G = crossprod(U, Vd)
  GWG = crossprod(Vd, K%*%Vd)
  Rmat=W%*%G%*%qr.solve(GWG, tol=qrtol)
  return(list(G=G, B=B, Rmat=Rmat, GWG=GWG))
} 
