## Julia Bischof
## 10-09-2015

clones.CDR3Length<-function(CDR3Length=NULL,functionality=NULL, junctionFr=NULL,abundance=c("relative","absolute"),...){
  if(length(CDR3Length)==0){
    stop("--> CDR3 length vector is missing")
  }
  if(length(abundance)!=1 || !(abundance %in% c("relative","absolute"))){
    abundance<-"relative"
  }
  
  out.list<-list()
  if(length(functionality)>0){
    functionality.new<-unlist(apply(data.frame(functionality),1,function(x){strsplit(x,split=" |,|[.]|;|[|]|_|-|/")[[1]]}))
    functionality.new<-functionality.new[which(nchar(functionality.new)>0)]
  }
  
  if(length(junctionFr)>0){
    junctionFr.new<-unlist(apply(data.frame(junctionFr),1,function(x){strsplit(x,split=" |,|[.]|;|[|]|_|/")[[1]]}))
    junctionFr.new<-junctionFr.new[which(nchar(junctionFr.new)>0)]
  }
  
    minCDR3<-min(as.numeric(CDR3Length),na.rm=T)
    maxCDR3<-max(as.numeric(CDR3Length),na.rm=T)
    # only CDR3 length
      tab.bar<-t(data.frame(apply(data.frame(seq(minCDR3,maxCDR3,1)),1,function(x){length(which(CDR3Length==x))})))
      colnames(tab.bar)<-paste("CDR3_length_",as.character(seq(minCDR3,maxCDR3,1)),sep="")
      if(abundance=="relative"){
        tab.bar<-tab.bar/length(CDR3Length)
      }
  if(length(functionality)==0 && length(junctionFr)==0){
    out.list<-data.frame(tab.bar,row.names=NULL)
  }else{
    out.list<-c(out.list,list(data.frame(tab.bar,row.names=NULL)))
    names(out.list)<-c(names(out.list)[which(names(out.list)!="")],"CDR3_length")
  }
  
  # CDR3 length vs. Functionality
    if(length(functionality)>0 && length(CDR3Length)==length(functionality)){
      CDR3Length.new<-vector()
      for(i in 1:length(CDR3Length)){
        CDR3Length.new<-c(CDR3Length.new,rep(CDR3Length[i],length(grep("prod|unknown|results",functionality[i]))))
      }
      tab.bar<-vector()
      for(i in c("^productive","^unproductive","unknown|No results")){
        tab.bar<-rbind(tab.bar,apply(data.frame(seq(minCDR3,maxCDR3,1)),1,function(x){length(intersect(grep(i,functionality.new),which(CDR3Length.new==x)))}))
      }
      colnames(tab.bar)<-paste("CDR3_length_",as.character(seq(minCDR3,maxCDR3,1)),sep="")
      rownames(tab.bar)<-c("productive","unproductive","unknown")
      if(abundance=="relative"){
        tab.bar<-t(data.frame(apply(tab.bar,1,function(x){x/colSums(tab.bar,na.rm=T)})))
        for(i in 1:nrow(tab.bar)){
          tab.bar[i,which(tab.bar[i,]=="NaN")]<-0
        }
      }
      out.list<-c(out.list,list(tab.bar))
      names(out.list)<-c(names(out.list)[which(names(out.list)!="")],"CDR3_length_vs_functionality")      
    }else if(length(functionality)>0 && length(CDR3Length)!=length(functionality)){
      stop("--> Lengths of CDR3 Length vector and functionality vector are different or missing")
    }
    
    # CDR3 length vs. JUNCTION_frame
    if(length(junctionFr)>0 && length(CDR3Length)==length(junctionFr)){
      CDR3Length.new<-vector()
      for(i in 1:length(CDR3Length)){
        CDR3Length.new<-c(CDR3Length.new,rep(CDR3Length[i],length(grep("frame|null",junctionFr[i]))))
      }
      tab.bar<-vector()
      for(i in c("in-frame","out-of-frame","null")){
        tab.bar<-rbind(tab.bar,apply(data.frame(seq(minCDR3,maxCDR3,1)),1,function(x){length(intersect(grep(i,junctionFr.new),which(CDR3Length.new==x)))}))
      }
      colnames(tab.bar)<-paste("CDR3_length_",as.character(seq(minCDR3,maxCDR3,1)),sep="")
      rownames(tab.bar)<-c("in-frame","out-of-frame","unknown")
      if(abundance=="relative"){
        tab.bar<-t(data.frame(apply(tab.bar,1,function(x){x/colSums(tab.bar,na.rm=T)})))
        for(i in 1:nrow(tab.bar)){
          tab.bar[i,which(tab.bar[i,]=="NaN")]<-0
        }
      }
      out.list<-c(out.list,list(tab.bar))
      names(out.list)<-c(names(out.list)[which(names(out.list)!="")],"CDR3_length_vs_junction_frame")
      
    }else if(length(junctionFr)>0 && length(CDR3Length)!=length(junctionFr)){
      stop("--> Lengths of 'CDR3 Length' vector and junction frame vector are different or missing")
    }
  return(out.list)
}


