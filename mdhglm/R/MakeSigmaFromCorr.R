MakeSigmaFromCorr <-
function(CorrMat=matrix(c(0,1,1,0),2,2),Correlations=c(0.5),Variances=c(1,1),indCorrIndex=10){
      
    TempCorrMat<-list(0)
    CorrMatOut<-list(0)
    TempCorrMat<-CorrMat
    for (j in 1:length(Correlations)){
        TempCorrMat[CorrMat==j]<-Correlations[j]
    }
    diag(TempCorrMat)<-1
    CorrMatOut<-TempCorrMat
    
    SigmaMat<-sqrt(Variances)*t(CorrMatOut*sqrt(Variances))
        
    SigmaTot<-SigmaMat%x%diag(indCorrIndex)
    invSigmaTot<-solve(SigmaMat)%x%diag(indCorrIndex)   

    return(list(SigmaTot=SigmaTot,invSigmaTot=invSigmaTot,SigmaMat=SigmaMat))
}
