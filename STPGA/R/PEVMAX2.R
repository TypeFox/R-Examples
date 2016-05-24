
PEVMAX2 <-
  function(Train,Test, P, lambda=1e-5){
    PTrain<-P[rownames(P)%in%Train,]
    PTest<-P[rownames(P)%in%Test,]
    PEV<-PTest%*%solve(crossprod(PTrain)+lambda*diag(ncol(P)),t(PTrain))
    PEVmean<-max(diag(tcrossprod(PEV)))
    return(PEVmean)
  }