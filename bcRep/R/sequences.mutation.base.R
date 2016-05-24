## Julia Bischof
## 2016-02-24

sequences.mutation.base<-function(mutationtab = NULL, summarytab = NULL, sequence=c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2"), nrCores=1){
  if(length(mutationtab)==0){
    stop("--> Mutationtab (7_V-REGION-mutation-and-AA-change-table(...).txt) is missing")
  }
  if(length(summarytab)==0){
    stop("--> Summarytab (1_Summary(...).txt) is missing")
  }
  if(length(sequence)!=1 && !(sequence %in% c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2"))){
    stop("--> Sequence is unknown or missing")
  }
  
  if(nrCores>4){
    nrCores<-4
  }
  cl<-makeCluster(nrCores)
  registerDoParallel(cl)
  
  
  mut.temp<-apply(data.frame(mutationtab[,grep(paste(sequence,"_REGION|",sequence,"_IMGT",sep=""), colnames(mutationtab))]),1,
                  function(x){strsplit(x, split=",|[\\|]|[(]")[[1]][grep("[a-z]",strsplit(x, split=",|[\\|]|[(]")[[1]])]})
  length.temp<-cbind(mutationtab$Sequence_ID, unlist(lapply(mut.temp,function(x){length(x)})))
  mut.temp<-unlist(mut.temp)[grep("[a-z]", unlist(mut.temp))]
  mut.list<-cbind(unlist(lapply(unlist(mut.temp),function(x){substr(x,1,1)})),
                  unlist(lapply(unlist(mut.temp),function(x){strsplit(x, split="[a-z]|[>]")[[1]][2]})))
  
  ids<-unlist(apply(length.temp,1, function(x){if(x[2]>0){rep(x[1],as.numeric(x[2]))}}))
  
  mut.list<-cbind(ids, mut.list)
  colnames(mut.list)<-c("Sequence_ID","from","position") #,"to")
  
  j<-NULL
  bases.ratio<-foreach(j=c("a","c","g","t")) %dopar% {
    cat(j, " ... ")
    subtab<-subset(mut.list, mut.list[,2]==j)
    subtab<-cbind(subtab,summarytab[match(subtab[,1], summarytab$Sequence_ID), "Sequence"])
    
    pos.m3<-apply(subtab,1,function(x){substr(x[4],as.numeric(x[3])-3,as.numeric(x[3])-3)})
    pos.m2<-apply(subtab,1,function(x){substr(x[4],as.numeric(x[3])-2,as.numeric(x[3])-2)})
    pos.m1<-apply(subtab,1,function(x){substr(x[4],as.numeric(x[3])-1,as.numeric(x[3])-1)})
    #pos.0<-apply(subtab,1,function(x){substr(x[4],as.numeric(x[3]),as.numeric(x[3]))})
    pos.p1<-apply(subtab,1,function(x){substr(x[4],as.numeric(x[3])+1,as.numeric(x[3])+1)})
    pos.p2<-apply(subtab,1,function(x){substr(x[4],as.numeric(x[3])+2,as.numeric(x[3])+2)})
    pos.p3<-apply(subtab,1,function(x){substr(x[4],as.numeric(x[3])+3,as.numeric(x[3])+3)})
    
    pos.m3<-pos.m3[-which(pos.m3=="")]
    pos.m2<-pos.m2[-which(pos.m2=="")]
    pos.m1<-pos.m1[-which(pos.m1=="")]
    #pos.0<-pos.0[-which(pos.0=="")]
    pos.p1<-pos.p1[-which(pos.p1=="")]
    pos.p2<-pos.p2[-which(pos.p2=="")]
    pos.p3<-pos.p3[-which(pos.p3=="")]
    
    
    pos.m3.ratio<-c(length(grep("a|A",pos.m3))/length(pos.m3),length(grep("c|C",pos.m3))/length(pos.m3),length(grep("g|G",pos.m3))/length(pos.m3),length(grep("t|T",pos.m3))/length(pos.m3))
    pos.m2.ratio<-c(length(grep("a|A",pos.m2))/length(pos.m2),length(grep("c|C",pos.m2))/length(pos.m2),length(grep("g|G",pos.m2))/length(pos.m2),length(grep("t|T",pos.m2))/length(pos.m2))
    pos.m1.ratio<-c(length(grep("a|A",pos.m1))/length(pos.m1),length(grep("c|C",pos.m1))/length(pos.m1),length(grep("g|G",pos.m1))/length(pos.m1),length(grep("t|T",pos.m1))/length(pos.m1))
    #pos.0.ratio<-c(length(grep("a|A",pos.0))/length(pos.0),length(grep("c|C",pos.0))/length(pos.0),length(grep("g|G",pos.0))/length(pos.0),length(grep("t|T",pos.0))/length(pos.0))
    pos.p1.ratio<-c(length(grep("a|A",pos.p1))/length(pos.p1),length(grep("c|C",pos.p1))/length(pos.p1),length(grep("g|G",pos.p1))/length(pos.p1),length(grep("t|T",pos.p1))/length(pos.p1))
    pos.p2.ratio<-c(length(grep("a|A",pos.p2))/length(pos.p2),length(grep("c|C",pos.p2))/length(pos.p2),length(grep("g|G",pos.p2))/length(pos.p2),length(grep("t|T",pos.p2))/length(pos.p2))
    pos.p3.ratio<-c(length(grep("a|A",pos.p3))/length(pos.p3),length(grep("c|C",pos.p3))/length(pos.p3),length(grep("g|G",pos.p3))/length(pos.p3),length(grep("t|T",pos.p3))/length(pos.p3))
    
    temp<-data.frame(pos.m3.ratio,pos.m2.ratio,pos.m1.ratio, rep(NA,4),pos.p1.ratio,pos.p2.ratio,pos.p3.ratio)
    colnames(temp)<-c("Pos-3","Pos-2","Pos-1","Pos0","Pos+1","Pos+2","Pos+3")
    rownames(temp)<-c(paste(j,"_to_a",sep=""),paste(j,"_to_c",sep=""),paste(j,"_to_g",sep=""),paste(j,"_to_t",sep=""))
    return(temp)
  }
  stopCluster(cl)
  
  bases.tab<-do.call(rbind.data.frame, bases.ratio)
  
  return(bases.tab)
}
