

PEVMAX <-
  function(Train,Test, P, lambda=1e-5){
    PTrain<-P[rownames(P)%in%Train,]
    PTest<-P[rownames(P)%in%Test,]
    PEVmax<-max(diag(PTest%*%solve(crossprod(PTrain)+lambda*diag(ncol(P)),t(PTest))))
    return(PEVmax)
  }
