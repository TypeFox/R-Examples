gcvssa <-
function(etas,yty,KtK,KtJ,Kty,Jty,JtJ,Qmats,ndpts,
         alpha,nknots,newlam,nbf,xnames,Jnames) {
  
  gammas <- exp(etas)
  gamvec <- NULL
  for(j in 1:length(Jnames)){
    xi <- strsplit(Jnames[j],":")
    xidx <- match(xi[[1]],xnames)
    gamvec <- c(gamvec,prod(gammas[xidx]))
  }
  gammat <- kronecker(gamvec,diag(nknots))
  jtj <- crossprod(gammat,JtJ)%*%gammat
  xty <- c(Kty,crossprod(gammat,Jty))
  ktj <- KtJ%*%gammat
  nqmat <- matrix(0,nbf+nknots,nbf+nknots)
  nqmat[(nbf+1):(nknots+nbf),(nbf+1):(nknots+nbf)] <- ndpts*newlam*(Qmats%*%gammat)
  xtx <- rbind(cbind(KtK,ktj),cbind(t(ktj),jtj))
  gcv <- tryCatch({
    chi <- pinvsm(xtx+nqmat)
    parta <- chi%*%xty
    gnum <- yty - 2*crossprod(xty,parta) + crossprod(parta,xtx%*%parta)
    ndpts*gnum/((ndpts-alpha*sum(diag(chi%*%xtx)))^2)
  }, error = function(e) yty)
  
}
