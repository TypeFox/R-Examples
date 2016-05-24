## Julia Bischof
## 10-09-2015

#library(doParallel)
#library(parallel)

sequences.geneComb<-function(family1=NULL,family2=NULL,level=c("subgroup", "gene", "allele"),abundance=c("relative", "absolute"),nrCores=1){
  if(length(family1)==0){
    stop("--> Vector of gene family 1 is missing")
  }
  if(length(family2)==0){
    stop("--> Vector of gene family 2 is missing")
  }
  if(length(level)!=1 || !(level %in% c("subgroup", "gene", "allele"))){
    stop("--> Gene level (subgroup, gene, allele) is missing")
  }
  if(length(abundance)!=1 || !(abundance %in% c("relative", "absolute"))){
    abundance<-"relative"
  }
   if(as.numeric(nrCores)>as.numeric(detectCores())){
    stop(paste("--> nrCores is higher than available number of cores (only ",as.numeric(detectCores())," cores available)",sep=""))
  }
  
  cl<-makeCluster(nrCores)
  registerDoParallel(cl)

    
  if(length(grep(" ",(family1[which(family1!="")[1]])))==0){
    familyname1<-substr(family1[which(family1!="")[1]],1,4)
  }else{
    familyname1<-substr(strsplit(family1[which(family1!="")[1]],split=" ")[[1]][2],1,4)
  }
    list1<-unique(unlist(apply(data.frame(family1),1,function(x){strsplit(x,split=" ")[[1]]})))
    list1<-unique(list1[grep(familyname1,list1)])
    if(level=="gene"){
      list1<-unique(unlist(apply(data.frame(list1),1,function(x){strsplit(x,split="[*]")[[1]][1]})))
    }else if(level=="subgroup"){
      list1<-unique(unlist(apply(data.frame(list1),1,function(x){strsplit(x,split="S|-|[*]")[[1]][1]})))
    }
  list1<-sort(list1)
  if(length(grep(" ",(family2[which(family2!="")[1]])))==0){
    familyname2<-substr(family2[which(family2!="")[1]],1,4)
  }else{
    familyname2<-substr(strsplit(family2[which(family2!="")[1]],split=" ")[[1]][2],1,4)
  }
  list2<-unique(unlist(apply(data.frame(family2),1,function(x){strsplit(x,split=" ")[[1]]})))
    list2<-unique(list2[grep(familyname2,list2)])
    if(level=="gene"){
      list2<-unique(unlist(apply(data.frame(list2),1,function(x){strsplit(x,split="[*]")[[1]][1]})))
    }else if(level=="subgroup"){
      list2<-unique(unlist(apply(data.frame(list2),1,function(x){strsplit(x,split="S|-|[*]")[[1]][1]})))
    }
    list2<-sort(list2)
  
    comb.names<-vector()
    comb.tab<-vector()
    list1.new<-c(list1,"")
    list2.new<-c(list2,"")
    i<-NULL
  if(level=="allele"){
    temp<-vector()
    comb.tab<-foreach(i=1:length(list1.new)) %dopar% {
      temp<-c(list1.new[i],unlist(apply(data.frame(list2.new),1,function(x){length(intersect(grep(gsub("[*]","_",x),gsub("[*]","_",family2),perl=T),grep(gsub("[*]","_",list1.new[i]), gsub("[*]","_",family1),perl=T)))})))
    }
  }else{
    temp<-vector()
    comb.tab<-foreach(i=1:length(list1.new)) %dopar% {
      temp<-c(list1.new[i],unlist(apply(data.frame(list2.new),1,function(x){length(intersect(grep(x,family2),grep(list1.new[i], family1)))})))
    }
  }
    
    comb.tab<-do.call(rbind.data.frame, comb.tab)
    newnames<-as.character(comb.tab[,1])
    comb.tab<-comb.tab[,-1]
    comb.tab<-apply(comb.tab,1,function(x){as.numeric(x)})
  
  colnames(comb.tab)<-c(as.character(newnames)[1:(length(newnames)-1)],"NA")
  rownames(comb.tab)<-c(as.character(list2),"NA")
    
    stopCluster(cl)
    
  if(abundance=='relative'){
      return(comb.tab/sum(comb.tab,na.rm=T))
    }
    return(comb.tab)
}


