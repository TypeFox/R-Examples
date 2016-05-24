ModelPCA<-function(iicc){
    
    nPCs<-iicc$parameterIdeal
    
    TFBS<-iicc$Transcriptionfactor
     Prob<-iicc$background
    missing<-iicc$missing
     NumericalMatrix<-numericalDNA(Prob)
    suma<-apply(TFBS,2,function(y){sum(y=="-")})
    threshold<-floor(nrow(TFBS)*missing/100)
    TFBS<-TFBS[, suma<=threshold]
    ncolTFBS<-ncol(TFBS)
    TFBSnum<-apply(TFBS,1,function(x){as.vector(t(NumericalMatrix[x,]))})
    TFBSnum<-t(TFBSnum)
       model<-pca(TFBSnum, nPcs=nPCs, method="svd", center=TRUE)
    JacksonPars<-JacksonParameters(nPCs,TFBSnum)
    parametersModel<-list(NumericalMatrix=NumericalMatrix,JacksonPars=JacksonPars, ncolTFBS=ncolTFBS)
    y<-list(model=model, parametersModel=parametersModel)
    return(y)
    
}
