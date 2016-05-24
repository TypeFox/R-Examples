gcvssp <-
function(etas,yty,KtK,KtJ,Kty,Jty,JtJ,Qmats,
         ndpts,alpha,nknots,newlam,nbf) {
  
  thmat <- kronecker(exp(etas),diag(nknots))
  jtj <- crossprod(thmat,JtJ)%*%thmat
  xty <- c(Kty,crossprod(thmat,Jty))
  ktj <- KtJ%*%thmat
  nqmat <- matrix(0,nknots+nbf,nknots+nbf)
  nqmat[(nbf+1):(nknots+nbf),(nbf+1):(nknots+nbf)] <- ndpts*newlam*(Qmats%*%thmat)
  xtx <- rbind(cbind(KtK,ktj),cbind(t(ktj),jtj))
  gcv <- tryCatch({
    chi <- pinvsm(xtx+nqmat)
    parta <- chi%*%xty
    gnum <- yty - 2*crossprod(xty,parta) + crossprod(parta,xtx%*%parta)
    ndpts*gnum/((ndpts-alpha*sum(diag(chi%*%xtx)))^2)
  }, error = function(e) yty)
  
}
