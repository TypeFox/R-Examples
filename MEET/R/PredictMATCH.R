PredictMATCH <-function(iicc){


    write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")

    Prob<-iicc$background
    Sequence <- iicc$DNA[[1]]
    ncolTFBS<-iicc$model$parametersModel$ncolTFBS
 
    Background<-matrix(,length(Sequence)-ncolTFBS+1,ncolTFBS)
    for(X in 1:(length(Sequence)-ncolTFBS+1)){
    Background[X,]<-Sequence[X:(X+ncolTFBS-1)]
  }
    
    CoreCut<-iicc$model$parametersModel$Corecut/100
    threshold<-iicc$threshold
    
      
      iicc$minim<-iicc$model$parametersModel$minim
      iicc$maxim<-iicc$model$parametersModel$maxim
      iicc$logodds<-iicc$model$parametersModel$logodds
      iicc$minimCore<-iicc$model$parametersModel$minim_core
      iicc$maximCore<-iicc$model$parametersModel$maxim_core
      iicc$core<-iicc$model$parametersModel$core
      iicc$index<-iicc$model$parametersModel$posCore
      CoreSeq<- Background[,iicc$index:(iicc$index+4)]
      
      Scores<-apply (Background,1, CalculScores,logodds=iicc$logodds)
        
      CoreScores<-apply(CoreSeq,1,CalculScores, logodds=iicc$core)
        
      Similarity<-sapply(Scores,CalculSimilarity, minim=iicc$minim, maxim=iicc$maxim)
      CoreSimilarity<-sapply(CoreScores,CalculSimilarity, minim=iicc$minimCore, maxim=iicc$maximCore)
      Scores<-cbind(Similarity, CoreSimilarity, deparse.level=0)
  
      DetectedFactors<-cbind(which(Scores[,1]>=threshold),Scores[Scores[,1]>=threshold,])
      DetectedFactors<-DetectedFactors[DetectedFactors[,2]>=CoreCut,]
      
      if(nrow(DetectedFactors)==0){
      output<-"No Binding Sites Found"
      }else{
      output<-lapply(c(1:nrow(DetectedFactors)),function(x){cbind(Sequence=paste(Sequence[DetectedFactors[x,1]:(DetectedFactors[x,1]+iicc$lenmotif)],sep="",collapse=""),pvalue=DetectedFactors[x,2], position=DetectedFactors[x,1])})
    }


}

