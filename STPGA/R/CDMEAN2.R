


CDMEAN2 <-
  function(Train,Test, P, lambda=1e-5){
    PTrain<-P[rownames(P)%in%Train,]
    PTest<-P[rownames(P)%in%Test,]
    CD<-PTest%*%solve(crossprod(PTrain)+lambda*diag(ncol(P)),t(PTrain))
    CDmean<-mean(diag(tcrossprod(CD))/diag(tcrossprod(PTest)))
    return(CDmean)
  }
