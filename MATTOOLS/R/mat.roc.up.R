 
  mat.roc.up<-  function (mattools.roclist,rocEvalSeq=seq(0,2,0.01)) 
  {
    tlen=length(mattools.roclist)-1
    for(i in 1:tlen){
      
      mattools.roclist[[i]][[3]]=mat.ROCcalc(truth=c(rep(0,length(mattools.roclist[[i]][2]$outzone)),rep(1,length(mattools.roclist[[i]][1]$zonezone))), data=c(mattools.roclist[[i]][2]$outzone,mattools.roclist[[i]][1]$zonezone), evalseq=rocEvalSeq)
      
    }
    mattools.roclist
  }