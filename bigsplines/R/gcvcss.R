gcvcss <-
function(eta,yty,xtx,xty,nqmat,ndpts,alpha) {
  
  gcv <- tryCatch({
    chi <- pinvsm(xtx+exp(eta)*nqmat)
    parta <- chi%*%xty
    gnum <- yty - 2*crossprod(xty,parta) + crossprod(parta,xtx%*%parta)
    ndpts*gnum/((ndpts-alpha*sum(diag(chi%*%xtx)))^2)
  }, error = function(e) yty)
  
}
