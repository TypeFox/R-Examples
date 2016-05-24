## Julia Bischof
## 10-09-2015

aaDistribution<-function(sequences=NULL,numberSeq=FALSE){
  if(length(sequences)==0 && length(grep("[A-Z]",sequences))==0){
    stop("--> sequences sequence vector is missing or empty")
  }
  allsequences<-unlist(apply(data.frame(sequences),1,function(x){strsplit(x,split=" |,|;|[.]|_|-|/")[[1]]}))
  allsequences<-allsequences[which(allsequences!="")]
  allsequenceslength<-unlist(apply(data.frame(allsequences),1,function(x){nchar(x)}))

  unisequenceslength<-sort(unique(allsequenceslength))
  uniAA<-c("F","L","I","M","V","S","P","T","A","Y","H","C","W","N","D","G","Q","R","K","E","*")
  
  AAdistr<-list()
  Numbersequences.tab<-vector()
  for(i in 1:length(unisequenceslength)){
    sequencesPosFreq.all<-vector()
    NumberSeq<-vector()
      samelength<-allsequences[which(allsequenceslength==unisequenceslength[i])]
      NumberSeq<-c(NumberSeq,length(which(allsequenceslength==unisequenceslength[i])))
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
    Numbersequences.tab<-rbind(Numbersequences.tab,NumberSeq)
  }  
  colnames(Numbersequences.tab)<-"Number_of_sequences"
  rownames(Numbersequences.tab)<-paste("sequence length = ",unisequenceslength,sep="")
  names(AAdistr)<-paste("sequence length = ",unisequenceslength,sep="")
  
  if(numberSeq==T){
    output<-list(AAdistr,Numbersequences.tab)
    names(output)<-c("Amino_acid_distribution","Number_of_sequences_per_length")
  }else{
    output<-list(AAdistr)
    names(output)<-"Amino_acid_distribution"
  }
  return(output)
}


