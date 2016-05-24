## Julia Bischof
## 15-10-2015

clones.filterSize<-function(clones.tab=NULL, column=NULL, number=NULL, propOfClones=NULL, propOfSequences=NULL, 
                            method=c("two.tailed","upper.tail","lower.tail")){
  if(length(clones.tab)==0){
    stop("--> Clones.tab is missing")
  }
  if(sum(c(length(number),length(propOfClones),length(propOfSequences)))>1 || sum(c(length(number),length(propOfClones),length(propOfSequences)))==0){
    stop("--> Please input one (!) filter criteria")
  }
  if(length(column)!=1){
    stop("--> Please insert column number of colname of column, which includes clone sizes (f.e. total_number_of_sequences)")
  }
  if(length(method)!=1 || !(method %in% c("two.tailed","upper.tail","lower.tail"))){
    stop("--> Method (two.tailed,upper.tail,lower.tail) is missing")
  }
  clones.tab<-clones.tab[order(as.numeric(clones.tab[,column]),decreasing = T),]
  
  if(length(number)>0){
    if(method=="two.tailed"){
      clones.up<-clones.tab[1:number,]
      clones.low<-clones.tab[(nrow(clones.tab)-number+1):nrow(clones.tab),]
    }else if(method=="upper.tail"){
      clones.up<-clones.tab[1:number,]
    }else if(method=="lower.tail"){
      clones.low<-clones.tab[(nrow(clones.tab)-number+1):nrow(clones.tab),]
    }
  }
  if(length(propOfClones)>0){
    tr<-round(propOfClones*nrow(clones.tab))
    if(method=="two.tailed"){
      clones.up<-clones.tab[1:tr,]
      clones.low<-clones.tab[(nrow(clones.tab)-tr+1):nrow(clones.tab),]
    }else if(method=="upper.tail"){
      clones.up<-clones.tab[1:tr,]
    }else if(method=="lower.tail"){
      clones.low<-clones.tab[(nrow(clones.tab)-tr+1):nrow(clones.tab),]
    }
  }
  if(length(propOfSequences)>0){
    tr<-as.numeric(sum(as.numeric(clones.tab[,column]),na.rm=T)*propOfSequences)
    if(method=="two.tailed"){
      clones.up<-clones.tab[which(clones.tab[,column]>tr),]
      clones.low<-clones.tab[which(clones.tab[,column]<=tr),]
    }else if(method=="upper.tail"){
      clones.up<-clones.tab[which(clones.tab[,column]>tr),]   
    }else if(method=="lower.tail"){
      clones.low<-clones.tab[which(clones.tab[,column]<=tr),]
    }
  }
  
  if(method=="two.tailed"){
    outlist<-list(clones.up, clones.low)
    names(outlist)<-c("upper.tail","lower.tail")
    return(outlist)
  }else if(method=="upper.tail"){
    return(clones.up)
  }else if(method=="lower.tail"){
    return(clones.low)
  }
}
