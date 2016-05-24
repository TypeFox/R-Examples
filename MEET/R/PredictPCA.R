PredictPCA <-function(iicc){


  require("pcaMethods")

  Prob<-iicc$background
  ncolTFBS<-iicc$model$parametersModel$ncolTFBS
  
  nPcs<-iicc$parametersIdeal
  threshold<-iicc$threshold
  Sequence <- iicc$DNA[[1]]
  
  missing<-iicc$missing
  NumericalMatrix<-iicc$model$parametersModel$NumericalMatrix
  model<-iicc$model$model
  parametersQ<-iicc$model$parametersModel$JacksonPars
  DNA<-sapply(X=c(1:(length(Sequence)-ncolTFBS+1)), FUN=function(X, ncolTFBS){Sequence[X:(X+ncolTFBS-1)]}, ncolTFBS=ncolTFBS)
  DNA<-t(DNA)
  DNAnum<-apply(DNA,1,function(x){as.vector(t(NumericalMatrix[x,]))})
  DNAnum<-t(DNAnum)
  matrix.residuals<-residuals(model,DNAnum, nPcs=model@nPcs)
  residus<-apply(matrix.residuals,1,function(vector){vector%*%vector})
  Score<-sapply(residus,QtoJackson, h0=parametersQ$h0, x1=parametersQ$x1,x2=parametersQ$x2,x3=parametersQ$x3 )
  
  DetectedFactors<-cbind(which((1-Score)<=threshold),Score[(1-Score)<=threshold])
 
   if(nrow(DetectedFactors)==0){
      output<-"No Binding Sites Found"
    }else{
      output<-lapply(c(1:nrow(DetectedFactors)),function(x){cbind(Sequence=paste(Sequence[DetectedFactors[x,1]:(DetectedFactors[x]+ncolTFBS)],sep="",collapse=""),pvalue=1-DetectedFactors[x,2], position=DetectedFactors[x,1], direction=iicc$direction)})
    }
  
   return(output)

}

