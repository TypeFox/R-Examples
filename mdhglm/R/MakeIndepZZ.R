MakeIndepZZ <-
function(ZZCorr,ZZmodel,Subject,SigmaMat=diag(2)){ 
    ZZCorrVec<-list(0)
    ncolZZCorr<-rep(0,length(ZZCorr))
    for (i in 1:length(ZZCorr)){
        ZZCorrVec[[i]]<-ZZCorr[[i]]%*%rep(1,ncol(ZZCorr[[i]]))
        ncolZZCorr[i] <-ncol(ZZCorr[[i]])
    }
               
    itchol<-t(chol(SigmaMat))  # This is actually cholesky decomposition instead of inverse, there was before inverse which was wrong
    CholeskyMatrices<-itchol
    maxModel<-as.numeric(max(names(table(ZZmodel))))
    ZZCorrUpd<-rep(list(0),maxModel)
    
    for (i in 1:maxModel) {         # loop over models 
         for (j in 1:length(ZZmodel)){ # loop over candidates
             if (ZZCorrUpd[[i]][1]==0 & length(ZZCorrUpd[[i]])==1) {
                 if (ZZmodel[j]==i) ZZCorrUpd[[i]]<-ZZCorrVec[[j]]
             }
             else {
                 if (ZZmodel[j]==i) ZZCorrUpd[[i]]<-cbind(ZZCorrUpd[[i]],ZZCorrVec[[j]])
             }
         }
    }
    for (i in 1:maxModel){
        if (i==1) ZZCorrUpdTot <- ZZCorrUpd[[1]]
        else ZZCorrUpdTot <- dbind(ZZCorrUpdTot, ZZCorrUpd[[i]])
   }
   
   ZZCorrUpdTotTemp <- ZZCorrUpdTot%*%itchol
   ZZOut<-rep(list(0),length(ZZCorr))
   for (i in 1:ncol(ZZCorrUpdTotTemp)) { 
        DiagDesign[[i]]<- model.matrix(~factor(Subject)-1)
        DiagDesign[[i]]<-DiagDesign[[i]]*ZZCorrUpdTotTemp[,i]
   }
    
    return(list(CholeskyMatrices=CholeskyMatrices, DiagDesign=DiagDesign))
    
}
