

#correct version but very close to the previous one. 
PEVMEAN2 <-
  function(Train,Test, P, lambda=1e-5){
    PTrain<-P[rownames(P)%in%Train,]
    PTest<-P[rownames(P)%in%Test,]
    PEV<-PTest%*%solve(crossprod(PTrain)+lambda*diag(ncol(P)),t(PTrain))
    PEVmean<-mean(diag(tcrossprod(PEV)))
    return(PEVmean)
  }

