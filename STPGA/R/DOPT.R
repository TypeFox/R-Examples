


DOPT <-
  function(Train,Test=NULL, P, lambda=1e-5){
    PTrain<-P[rownames(P)%in%Train,]
    D<--det(crossprod(PTrain)+lambda*diag(ncol(P)))
    return(D)
  }