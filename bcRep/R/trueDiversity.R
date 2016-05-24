## Julia Bischof
## 10-09-2015

#library(vegan)

trueDiversity<-function(sequences=NULL,aaDistribution.tab=NULL,order=c(0, 1, 2)){
  
  if(length(sequences)==0 && length(aaDistribution.tab)==0){
    stop("--> Need either sequences or aaDistribution.tab as input")
  }
  
  if(length(order)!=1 || !(order %in% c(0, 1, 2))){
    stop("--> Diversity order is missing (0, 1 or 2)")
  }
  if(length(aaDistribution.tab)>0){
    AAdistr<-aaDistribution.tab$Amino_acid_distribution
    unisequenceslength<-as.numeric(apply(data.frame(names(aaDistribution.tab)),1,function(x){strsplit(x,split=" = ")[[1]][2]}))
  }else if(length(sequences)>0){
    allsequences<-unlist(apply(data.frame(sequences),1,function(x){strsplit(x,split=" |,|;|[.]|_|-")[[1]]}))
    allsequences<-allsequences[which(allsequences!="")]
    allsequenceslength<-unlist(apply(data.frame(allsequences),1,function(x){nchar(x)}))
    
    unisequenceslength<-sort(unique(allsequenceslength))
    uniAA<-c("F","L","I","M","V","S","P","T","A","Y","H","C","W","N","D","G","Q","R","K","E","*")
    
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
  }else if(length(sequences)>0){
    stop("--> Need either sequences sequence vector or Amino acid distribution tab as input")
  }
  
  truDiv<-list()                         
  for(i in 1:length(AAdistr)){
    if(length(which(is.na(AAdistr[[i]])))!=length(AAdistr[[i]])){
      AAdistr[[i]]<-apply(AAdistr[[i]],1,function(x){as.numeric(x)})
      if(order==0){
        td<-t(data.frame(specnumber(AAdistr[[i]])))
      }else if(order==1){
        td<-t(data.frame(exp(diversity(AAdistr[[i]],index="shannon",MARGIN=1,base=exp(1)))))
      }else if(order==2){
        td<-t(data.frame(diversity(AAdistr[[i]],index="invsimpson",MARGIN=1)))   
      }
    }else{
      td<-t(data.frame(rep(NA,length(AAdistr[[i]]))))
    }
    td<-as.data.frame(td,row.names="")
    colnames(td)<-paste("Position",seq(1,length(td),1),sep="")
    truDiv<-c(truDiv,list(td))
  }
  names(truDiv)<-names(AAdistr)
  truDiv<-c(list(order,truDiv))
  names(truDiv)<-c("True_diversity_order","True_diversity")
  return(truDiv)
}

