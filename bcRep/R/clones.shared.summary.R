## Julia Bischof
## 10-09-2015

clones.shared.summary<-function(shared.tab=NULL, clones.tab=NULL){
  if(length(shared.tab)==0){
    stop("--> Output of clones.shared() needed as input")
  }
  if(is.data.frame(shared.tab)==T){
    shind<-as.vector(shared.tab[,2])
  }else{
    shind<-as.vector(shared.tab)
  }
  
  unishind<-unique(shind)
  nrshclones<-vector() # shared clones
  for(i in 1:length(unishind)){
    nrshclones<-c(nrshclones,length(which(shind==unishind[i])))
  }
  nrshclones<-cbind(as.character(unishind),nrshclones)
  
  if(length(clones.tab)>0){
    if(is.data.frame(clones.tab)==T){
      clones.tab<-as.vector(clones.tab[,1])
    }else{
      clones.tab<-as.vector(clones.tab)
    }
    ind<-unique(clones.tab)
    
    clones.nr.new<-vector() # number clones per individual minus shared clones
    for(i in 1:length(ind)){
      clones.nr.new<-c(clones.nr.new,length(which(clones.tab==ind[i]))-sum(as.numeric(nrshclones[grep(ind[i],nrshclones[,1]),2])))
    }
    nrshclones<-rbind(cbind(paste("only in ",as.character(ind),sep=""),clones.nr.new),nrshclones)
    colnames(nrshclones)<-c("group",'number_clones')
  }
  return(nrshclones)
}


