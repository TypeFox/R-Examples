## Julia Bischof
## 10-15-2015

#library(doParallel)
#library(parallel)

compare.aaDistribution<-function(sequence.list=NULL, names=NULL, numberSeq=FALSE, nrCores=1){
  if(length(sequence.list)<2 && is.list(sequence.list)==F){
    stop("--> Need a list of at least 2 vectors")
  }
  
  if(as.numeric(nrCores)>as.numeric(detectCores())){
    stop(paste("--> nrCores is higher than available number of cores (only ",as.numeric(detectCores())," cores available)",sep=""))
  }
  
  cl<-makeCluster(nrCores)
  registerDoParallel(cl)

  
  uniAA<-c("F","L","I","M","V","S","P","T","A","Y","H","C","W","N","D","G","Q","R","K","E","*")
  
  AAdistr<-list()
  j<-NULL
  AAdistr.tab<-foreach(j=1:length(sequence.list)) %dopar% {
    outlist<-vector()
    allsequences<-unlist(apply(data.frame(sequence.list[[j]]),1,function(x){strsplit(x,split=" |,|;|[.]|_|-|/")[[1]]}))
    allsequences<-allsequences[which(allsequences!="")]
    allsequenceslength<-unlist(apply(data.frame(allsequences),1,function(x){nchar(x)}))
    
    unisequenceslength<-sort(unique(allsequenceslength))
    
    AAdistr<-list()
    for(i in 1:length(unisequenceslength)){
      sequencesPosFreq.all<-vector()
      samelength<-allsequences[which(allsequenceslength==unisequenceslength[i])]
      if(length(samelength)>0){
        if(unisequenceslength[i]==1){
          sequencesPos<-data.frame(apply(data.frame(samelength),1,function(x){strsplit(x,split="")[[1]]}))
        }else{
          sequencesPos<-t(apply(data.frame(samelength),1,function(x){strsplit(x,split="")[[1]]}))
        }
        sequencesPosFreq<-vector()
        for(k in 1:ncol(sequencesPos)){
          sequencesPosFreq<-rbind(sequencesPosFreq,apply(data.frame(uniAA),1,function(x){if(length(which(sequencesPos[,k]==x))>0){length(which(sequencesPos[,k]==x))}else{0}}))    
        }
        sequencesPosFreq<-t(sequencesPosFreq)
        sequencesPosFreq.all<-rbind(sequencesPosFreq.all,sequencesPosFreq)
      }else{
        sequencesPosFreq.all<-rbind(sequencesPosFreq.all,rep(NA,unisequenceslength[i]))
      }
      sequencesPosFreq.all<-sequencesPosFreq.all/colSums(sequencesPosFreq.all,na.rm=T)
      colnames(sequencesPosFreq.all)<-paste("Position",seq(1,ncol(sequencesPosFreq.all),1),sep="")
      rownames(sequencesPosFreq.all)<-uniAA
      AAdistr<-c(AAdistr,list(sequencesPosFreq.all))
    }  
    names(AAdistr)<-paste("sequence length = ",unisequenceslength,sep="")   
    s<-j
    AAdistr<-c(AAdistr, s)
    outlist<-AAdistr
  }
  
  index<-vector()
  for(i in 1:length(AAdistr.tab)){
    index<-c(index,AAdistr.tab[[i]][length(AAdistr.tab[[i]])])
    AAdistr.tab[[i]]<-AAdistr.tab[[i]][-length(AAdistr.tab[[i]])]
  }
  if(length(names)==length(sequence.list)){
    names(AAdistr.tab)<-names[sort(as.numeric(index))]
  }else{
    newnames<-paste("Sample",1:length(sequence.list),sep="")
    names(AAdistr.tab)<-newnames[sort(as.numeric(index))]
  }
  
  stopCluster(cl)
  
  if(numberSeq==T){
  
   cl<-makeCluster(nrCores)
   registerDoParallel(cl)

    Numbersequences.tab<-vector()
    NumberSeq<-vector()
    j<-NULL
    Numbersequences.tab<-foreach(j=1:length(sequence.list)) %dopar% {
      outtab<-vector()
      allsequences<-unlist(apply(data.frame(sequence.list[[j]]),1,function(x){strsplit(x,split=" |,|;|[.]|_|-|/")[[1]]}))
      allsequences<-allsequences[which(allsequences!="")]
      allsequenceslength<-unlist(apply(data.frame(allsequences),1,function(x){nchar(x)}))
      
      unisequenceslength<-sort(unique(allsequenceslength))
      
      for(i in 1:length(unisequenceslength)){
        NumberSeq<-cbind(NumberSeq,length(which(allsequenceslength==unisequenceslength[i])))
      }
      NumberSeq<-data.frame(NumberSeq,row.names=NULL)
      colnames(NumberSeq)<-paste("sequence length = ",unisequenceslength,sep="")
      s<-j
      outtab<-c(NumberSeq,s)
    }
    
    index<-vector()
    for(i in 1:length(Numbersequences.tab)){
      index<-c(index,Numbersequences.tab[[i]][length(Numbersequences.tab[[i]])])
      Numbersequences.tab[[i]]<-Numbersequences.tab[[i]][-length(Numbersequences.tab[[i]])]
    }
    
    if(length(names)==length(sequence.list)){
      names(Numbersequences.tab)<-names[sort(as.numeric(index))]
    }else{
      newnames<-paste("Sample",1:length(sequence.list),sep="")
      names(Numbersequences.tab)<-newnames[sort(as.numeric(index))]
    }
    
    stopCluster(cl)
    
    output<-list(AAdistr.tab,Numbersequences.tab)
    names(output)<-c("Amino_acid_distribution","Number_of_sequences_per_length")
  }else{
    output<-list(AAdistr.tab)
    names(output)<-"Amino_acid_distribution"
  }
  return(output)
  
}
