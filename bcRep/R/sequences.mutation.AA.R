## Julia Bischof
## 2016-02-24

sequences.mutation.AA<-function(mutationtab = NULL, sequence=c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2")){
  
  if(length(mutationtab)==0){
    stop("--> Mutationtab (7_V-REGION-mutation-and-AA-change-table(...).txt) is missing")
  }
  if(length(sequence)!=1 && !(sequence %in% c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2"))){
    stop("--> Sequence is unknown or missing")
  }
  
  AArange<-c("I", "V", "L", "F", "C", "M", "A", "W", 
             "G", "T", "S", "Y", "P", "H",
             "N", "D", "Q", "E", "K", "R")
  
  mut.list<-unlist(apply(data.frame(mutationtab[,grep(paste(sequence,"_REGION|",sequence,"_IMGT",sep=""), colnames(mutationtab))]),1,
                         function(x){strsplit(x, split=",|[\\|]|[(]")}))
  mut.list<-mut.list[grep("[A-Z]", mut.list)]
  mut.list<-paste(lapply(mut.list,function(x){substr(x,1,1)}),lapply(mut.list, function(x){strsplit(x, split="[>]")[[1]][2]}),sep="_")
  mut.uni<-unique(mut.list)
  temp<-unlist(lapply(mut.uni,function(x){length(which(mut.list==x))/length(mut.list)}))
  prop.tab<-matrix(0, ncol=length(AArange), nrow=length(AArange))
  for(i in 1:length(AArange)){
    for(j in 1:length(AArange)){
      if(length(which(mut.uni==paste(AArange[i],AArange[j],sep="_")))>0){
        prop.tab[i,j]<-temp[which(mut.uni==paste(AArange[i],AArange[j],sep="_"))]  
      }else{
        prop.tab[i,j]<-0
      }
    }
  }
  colnames(prop.tab)<-AArange # to
  rownames(prop.tab)<-AArange # from
  
  return(prop.tab)
}