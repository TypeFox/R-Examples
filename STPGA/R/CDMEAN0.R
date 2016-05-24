

CDMEAN0 <-
  function(Train,Test, P, lambda=1e-5){
    PTrain<-P[rownames(P)%in%Train,]
    CDmean<-mean(diag(PTrain%*%solve(crossprod(PTrain)+lambda*diag(ncol(P)),t(PTrain)))/diag(tcrossprod(PTrain)))
    return(CDmean)
  }
