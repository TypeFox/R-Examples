## Julia Bischof
## 10-15-2015

#library(vegan)
#library(doParallel)


compare.trueDiversity<-function(sequence.list=NULL,comp.aaDistribution.tab=NULL,order=c(0, 1, 2), 
                               names=NULL, nrCores=1){
  if(length(sequence.list)==0 && length(comp.aaDistribution.tab)==0){
    stop("--> Need either sequences or comp.aaDistribution.tab as input")
  }
  if(length(sequence.list)<2 && length(comp.aaDistribution.tab$Amino_acid_distribution)<2){
    stop("--> Need a list of at least 2 vectors")
  }
  
  if(length(order)!=1 || !(order %in% c(0, 1, 2))){
    stop("--> Diversity order is missing (0, 1 or 2)")
  }
  if(length(comp.aaDistribution.tab)>0){
    AAdistr<-comp.aaDistribution.tab$Amino_acid_distribution  
  }else if(length(sequence.list)>0 && length(comp.aaDistribution.tab)==0){
    temp<-sequence.list
    names.temp<-names
    nrCores.temp<-nrCores
    AAdistr<-compare.aaDistribution(sequence.list=temp, numberSeq = F, names = names.temp, nrCores=nrCores.temp)
    AAdistr<-AAdistr$Amino_acid_distribution
  }
  
  unisequenceslength<-vector()
  for(i in 1:length(AAdistr)){
    unisequenceslength<-unique(c(unisequenceslength,as.numeric(apply(data.frame(names(AAdistr[[i]])),1,function(x){strsplit(x,split=" = ")[[1]][2]}))))
  }
  
  truDiv.list<-list(order)
  for(j in 1:length(AAdistr)){
    truDiv<-list()                         
    for(i in 1:length(AAdistr[[j]])){
      if(length(which(is.na(AAdistr[[j]][[i]])))!=length(AAdistr[[j]][[i]])){
        AAdistr[[j]][[i]]<-apply(AAdistr[[j]][[i]],1,function(x){as.numeric(x)})
        if(order==0){
          td<-t(data.frame(specnumber(AAdistr[[j]][[i]])))
        }else if(order==1){
          td<-t(data.frame(exp(diversity(AAdistr[[j]][[i]],index="shannon",MARGIN=1,base=exp(1)))))
        }else if(order==2){
          td<-t(data.frame(diversity(AAdistr[[j]][[i]],index="invsimpson",MARGIN=1)))   
        }
      }else{
        td<-t(data.frame(rep(NA,length(AAdistr[[j]][[i]]))))
      }
      td<-as.data.frame(td,row.names="")
      colnames(td)<-paste("Position",seq(1,length(td),1),sep="")
      truDiv<-c(truDiv,list(td))
    }
    names(truDiv)<-names(AAdistr[[j]])
    #truDiv<-c(list(order,truDiv))
    #names(truDiv)<-c("True_diversity_order","True_diversity")
    truDiv.list<-c(truDiv.list,list(truDiv))
  }
  names(truDiv.list)<-c("True_diversity_order",if(length(names)==length(AAdistr)){names}else if(length(comp.aaDistribution.tab$Amino_acid_distribution)>0){names(AAdistr)}else{paste("Sample",1:length(AAdistr),sep="")})
  return(truDiv.list)  
}
