#####Summarize the accuracy of nodes and networking
SummaryAccuracy=function(Trace,No.Sets,Type.Accuracy=c("Node","Sub-network"),True.Node, 
                         TruePositive.Net,FalsePositive.Net,
                         Type.Net.Accuracy=c("Marginal","Sample"),Tolerance)
{
  
  results=list()
  if(Type.Accuracy=="Node")
  {
    indexone=which(True.Node==1)
    indexzero=which(True.Node==0)
    ###Checking True positive
    finalratioTsum=c()
    for (num in indexone)
    {
      finalratioT=TPR(num,Trace)
      finalratioTsum=c(finalratioTsum,finalratioT)
    }
    
    ###Checking False positive
    finalratioFsum=c()
    for (num in indexzero)
    {
      finalratioF=TPR(num,Trace)
      finalratioFsum=c(finalratioFsum,finalratioF)
    }
    results$TPR=finalratioTsum
    results$TPR.average=mean(finalratioTsum)
    results$FPR=finalratioFsum
    results$FPR.average=mean(finalratioFsum)
    results$FDR=(1-results$TPR.average)*length(indexone)/(length(indexone)+length(indexzero))+results$FPR.average*length(indexzero)/(length(indexone)+length(indexzero))
    show(results$TPR.average)
    show(results$FPR.average)
    show(results$FDR.average)
  }else if(Type.Accuracy=="Sub-network")
  {
    if (Type.Net.Accuracy=="Marginal")
    {
      x=0
      y=0
      for (i in 0:(No.Sets))
      {
        cat("Data Set: ",i,"\n")
        flush.console()
        
        j=c((i*nrow(Trace)/No.Sets+1):(i*nrow(Trace)/No.Sets+nrow(Trace)/No.Sets))
        new1=Trace[j,]
        ztall=sapply(TruePositive.Net,function(kk) return(mean(new1[,kk])))
        zfall=sapply(FalsePositive.Net,function(kk) return(mean(new1[,kk])))
        if(length(which(ztall>0.5))>=length(TruePositive.Net)*Tolerance & length(which(zfall>0.5))<=length(FalsePositive.Net)*(1-Tolerance)){x=x+1}
        if(length(which(ztall>0.5))>=length(TruePositive.Net)*Tolerance & length(which(zfall>0.5))>=length(FalsePositive.Net)*(1-Tolerance)){y=y+1}
        
      }
      results$TPR.average=x/No.Sets
      results$FPR.average=y/No.Sets
      results$FDR.average=y/(x+y)
      show(results)
    }else if (Type.Net.Accuracy=="Sample")
    {
      x=0
      y=0
      for (i in 1:nrow(Trace))
      {
        if (i%%(nrow(Trace)/No.Sets)==0){
          cat("Data Set: ",i/(nrow(Trace)/No.Sets),"\n")
          flush.console()}
        if (length(which(Trace[i,TruePositive.Net]>=1))>=length(TruePositive.Net)*Tolerance & length(which(Trace[i,FalsePositive.Net]>=1))<=length(FalsePositive.Net)*(1-Tolerance)){x=x+1}
        if (length(which(Trace[i,TruePositive.Net]>=1))>=length(TruePositive.Net)*Tolerance & length(which(Trace[i,FalsePositive.Net]>=1))>=length(FalsePositive.Net)*(1-Tolerance)){y=y+1}
      }
      results$TPR.average=x/nrow(Trace)
      results$FPR.average=y/nrow(Trace)
      results$FDR.average=y/(x+y)
      show(results)
    }
    
  }
  invisible(results)
}
