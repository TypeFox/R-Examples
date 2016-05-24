lamloop <-
function(lambdas,thvec,Kty,Jty,KtK,KtJ,JtJ,
         Qmats,nknots,ndpts,alpha,yty,nbf) {
  lamgcv <- vector("numeric",length(lambdas))
  nqmat <- matrix(0,nbf+nknots,nbf+nknots)
  thmat <- kronecker(thvec,diag(nknots))
  jtj <- crossprod(thmat,JtJ)%*%thmat
  xty <- c(Kty,crossprod(thmat,Jty))
  ktj <- KtJ%*%thmat
  nqmat[(nbf+1):(nknots+nbf),(nbf+1):(nknots+nbf)] <- ndpts*(Qmats%*%thmat)
  xtx <- rbind(cbind(KtK,ktj),cbind(t(ktj),jtj))
  for(jj in 1:length(lambdas)){
    lamgcv[jj] <- tryCatch({
      chi <- pinvsm(xtx+lambdas[jj]*nqmat)
      parta <- chi%*%xty
      gnum <- yty - 2*crossprod(xty,parta) + crossprod(parta,xtx%*%parta)
      ndpts*gnum/((ndpts-alpha*sum(diag(chi%*%xtx)))^2)
    }, error = function(e) yty)
  }
  newlam <- lambdas[which.min(lamgcv)]
}
