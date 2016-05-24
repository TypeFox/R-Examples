lamcoef <-
function(lambdas,thvec,Kty,Jty,KtK,KtJ,JtJ,
         Qmats,nknots,ndpts,alpha,yty,nbf){
  lamgcv <- vector("numeric",length(lambdas))
  ssval <- vector("numeric",length(lambdas))
  trval <- vector("numeric",length(lambdas))
  isqrt <- vector("list",length(lambdas))
  fcoefs <- matrix(0,nknots+nbf,length(lambdas))
  nqmat <- matrix(0,nbf+nknots,nbf+nknots)
  thmat <- kronecker(thvec,diag(nknots))
  jtj <- crossprod(thmat,JtJ)%*%thmat
  xty <- c(Kty,crossprod(thmat,Jty))
  ktj <- KtJ%*%thmat
  nqmat[(nbf+1):(nknots+nbf),(nbf+1):(nknots+nbf)] <- ndpts*(Qmats%*%thmat)
  xtx <- rbind(cbind(KtK,ktj),cbind(t(ktj),jtj))
  for(jj in 1:length(lambdas)){
    lamgcv[jj] <- tryCatch({
      ceig <- eigen(xtx+lambdas[jj]*nqmat,symmetric=TRUE)
      nze <- sum(ceig$val>ceig$val[1]*.Machine$double.eps)
      isqrt[[jj]] <- ceig$vec[,1:nze]%*%diag(ceig$val[1:nze]^-0.5)
      chi <- tcrossprod(isqrt[[jj]])
      fcoefs[,jj] <- chi%*%xty
      trval[jj] <- sum(diag(chi%*%xtx))
      ssval[jj] <- yty - 2*crossprod(xty,fcoefs[,jj]) + crossprod(fcoefs[,jj],xtx%*%fcoefs[,jj])
      ndpts*ssval[jj]/((ndpts-alpha*trval[jj])^2)
    }, error = function(e) yty)
  }
  opti <- which.min(lamgcv)
  fxinfo <- list(c(fcoefs[,opti],ssval[opti],trval[opti],opti),isqrt[[opti]])
}
