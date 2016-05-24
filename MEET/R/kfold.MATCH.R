kfold.MATCH <-function(iicc, Seqin){
    require("MEET")
    write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")
  
  position<-iicc$position
  vectorP<-iicc$vector
  
  Prob<-iicc$background

	  
	  x<-read.fasta(file=Seqin)
	  y<-x[-iicc$outsequence]
	      
       
    DNA <- iicc$DNAreal
  
  ROCcurve<-vector(mode="list", length=length(vectorP))

  
  for(i in 1:length(vectorP)){
    CoreCut<-vectorP[i]
    
    matriu_roc <- matrix(,,2)
    for(m in 1:length(y)){

	
	    z<-y[-m]
	    write.fasta(sequences=z,names=c(1:length(z)),file.out="background.fa",open="w") 
	    Factortrans<-switch(iicc$alignment, "CLUSTALW"=align.clustalw(filein="SetTF.fa", fileout="Sq.fa", call=iicc$call.clustalw), "MUSCLE"=align.muscle(filein="SetTF.fa", fileout="Sq.fa", gapopen=iicc$gapopen, maxiters=iicc$maxiters,gapextend=iicc$gapextend, call=iicc$call.muscle),"MEME"=align.MEME(filein="SetTF.fa",fileout="Sq.fa",iicc),"NONE"=Read.aligned("background.fa"), stop("Alignment method not included"))
	    
	
	window<-ncol(Factortrans) 
	suma<-apply(Factortrans,2,function(y){sum(y=="-")})
	training.set<-Factortrans<-Factortrans[, suma==0]
	iicc$Transcriptionfactor<-Factortrans
	
	 Background<-sapply(X=c(1:(length(DNA)-ncol(Factortrans)+1)), FUN=function(X, Factortrans){DNA[X:(X+ncol(Factortrans)-1)]}, Factortrans=Factortrans)
	Background<-t(Background)
	
	logodds <- CalculInformation(training.set, Prob=Prob)
	vector_informacio <- as.vector(logodds[,1])
	inf <- vector(mode='logical', length=(length(vector_informacio)-4))
	for (j in 1:(length(vector_informacio)-4)){
	    inf[j] <- vector_informacio[j]+vector_informacio[j+1]+vector_informacio[j+2]+vector_informacio[j+3]+vector_informacio[j+4]
	  }
	index <- which.max(inf)
	core <- logodds[(index:(index+4)),2:ncol(logodds)]
	logodds <- logodds[,2:dim(logodds)[2]]
	 minim<-0
      maxim<-0
      minim_core<-0
      maxim_core<-0
      for(j in 1:dim(logodds)[1]){
	  minim <- min(logodds[j,])+minim
	  maxim <- max(logodds[j,])+maxim
	  }
      for (j in 1:dim(core)[1]){
	    minim_core <- minim_core+ min(core[j,])
	    maxim_core <- maxim_core+ max(core[j,])
      }
      
      iicc$minim<-minim
     
      iicc$maxim<-maxim
     
      iicc$logodds<-logodds
     
      iicc$minimCore<-minim_core
     
      iicc$maximCore<-maxim_core
     
      iicc$core<-core
      iicc$index<-index
      CoreSeq<- Background[,index:(index+4)]
      
      Scores<-apply (Background,1, CalculScores,logodds=iicc$logodds)
      CoreScores<-apply(CoreSeq,1,CalculScores, logodds=iicc$core)
     
      Similarity<-sapply(Scores,CalculSimilarity, minim=iicc$minim, maxim=iicc$maxim)
      
      CoreSimilarity<-sapply(CoreScores,CalculSimilarity, minim=iicc$minimCore, maxim=iicc$maximCore)

      matriuROC<-rbind(Similarity, CoreSimilarity)
      matriuROC<-as.matrix(matriuROC)
      matriuROC<-t(matriuROC)
     
      matriu_roc<-rbind(matriu_roc, matriuROC)
    }
      matriu_roc<-matriu_roc[2:nrow(matriu_roc),]
    
     indexs <- which(matriu_roc[,2]<=CoreCut)
      print(summary(indexs))
     matriu_roc[indexs,1]<-0
       
    
    
     ROCcurve[[i]]<-matriu_roc[,1]
     

  }

  ROCcurve

}

