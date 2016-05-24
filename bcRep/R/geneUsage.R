## Julia Bischof
## 10-09-2015

geneUsage<-function(genes=NULL,level=c("subgroup", "gene", "allele"),functionality=NULL,junctionFr=NULL,
                           abundance=c("relative","absolute"),...){
  if(length(genes)==0){
    stop("--> Gene usage vector is missing")
  }
  if(length(level)!=1 || !(level %in% c("subgroup", "gene", "allele"))){
    stop("--> Level is missing (subgroup, gene or allele)")
  }
  if(level=="allele" && length(grep("[*]",genes))==0){
    stop("--> No allele information available in input vector")
  }
  if(length(abundance)!=1 || !(abundance %in% c("relative","absolute"))){
    abundance<-"relative"
  }
  
  out.list<-list()
  
  if(length(genes)>0){
    if(length(grep(" ",(genes[1])))==0){
      family<-substr(genes[1],1,4)
    }else{
      family<-substr(strsplit(genes[1],split=" ")[[1]][2],1,4)
    }
    genelist<-unlist(apply(data.frame(genes),1,function(x){strsplit(x,split=" |,|[.]|;|[|]|_|/")[[1]]}))
    genelist<-sort(genelist[grep(family,genelist)])
    if(level=="subgroup"){
      genelist<-sort(unlist(apply(data.frame(genelist),1,function(x){strsplit(x,split="S|-|[*]")[[1]][1]})))
    }else if(level=="gene"){
      genelist<-sort(unlist(apply(data.frame(genelist),1,function(x){strsplit(x,split="[*]")[[1]][1]})))
    }
    
    # only gene usage
    tab.bar<-t(data.frame(apply(data.frame(unique(genelist)),1,function(x){length(which(genelist==x))})))
    colnames(tab.bar)<-unique(genelist)
    if(abundance=="relative"){
      tab.bar<-tab.bar/length(genelist)
    }
    tab.bar<-data.frame(tab.bar,row.names=NULL)
    colnames(tab.bar)<-unique(genelist)
    if(length(functionality)==0 && length(junctionFr)==0){
      out.list<-tab.bar
    }
    out.list<-c(out.list,list(data.frame(tab.bar,row.names=NULL, check.names=F)))
    names(out.list)<-c(names(out.list)[which(names(out.list)!="")],"gene_usage")
    
    # gene usage vs. functionality
    if(length(functionality)>0){
      functionality.new<-unlist(apply(data.frame(functionality),1,function(x){strsplit(x,split=",|[.]|;|[|]|_|/")[[1]]}))
      functionality.new<-functionality.new[which(nchar(functionality.new)>0)]
      
      genelist.new<-vector()
      for(i in 1:length(genelist)){
        genelist.new<-c(genelist.new,rep(genes[i],length(grep("prod|unknown|results",functionality[i]))))
      }      
      
      tab.bar<-vector()
      for(i in c("^productive","^unproductive","unknown|No results")){
        tab.bar<-rbind(tab.bar,apply(data.frame(unique(genelist)),1,function(x){length(intersect(grep(i,functionality.new),grep(paste(x,"-|",x,"[*]|",x,"S|",gsub("[*]","_",x),sep=""),gsub("[*]","_",genelist.new))))}))
      }
      colnames(tab.bar)<-unique(genelist)
      rownames(tab.bar)<-c("productive","unproductive","unknown")
      if(abundance=="relative"){
        tab.bar<-t(data.frame(apply(tab.bar,1,function(x){x/colSums(tab.bar,na.rm=T)})))
        for(i in 1:nrow(tab.bar)){
          tab.bar[i,which(tab.bar[i,]=="NaN")]<-0
        }
      }
      out.list<-c(out.list,list(data.frame(tab.bar,check.names=F)))
      names(out.list)<-c(names(out.list)[which(names(out.list)!="")],"gene_usage_vs_functionality")      
    }
    
    # gene usage vs. JUNCTION_frame
    if(length(junctionFr)>0){
      junctionFr.new<-unlist(apply(data.frame(junctionFr),1,function(x){strsplit(x,split=" |,|[.]|;|[|]|_|/")[[1]]}))
      junctionFr.new<-junctionFr.new[which(nchar(junctionFr.new)>0)]
      
      genelist.new<-vector()
      for(i in 1:length(genelist)){
        genelist.new<-c(genelist.new,rep(genes[i],length(grep("frame|null",junctionFr[i]))))
      }        
      tab.bar<-vector()
      for(i in c("in-frame","out-of-frame","null")){
        tab.bar<-rbind(tab.bar,apply(data.frame(unique(genelist)),1,function(x){length(intersect(grep(i,junctionFr.new),grep(paste(x,"-|",x,"[*]|",x,"S|",gsub("[*]","_",x),sep=""),gsub("[*]","_",genelist.new))))}))
      }
      colnames(tab.bar)<-unique(genelist)
      rownames(tab.bar)<-c("in-frame","out-of-frame","unknown")
      if(abundance=="relative"){
        tab.bar<-t(data.frame(apply(tab.bar,1,function(x){x/colSums(tab.bar,na.rm=T)})))
        for(i in 1:nrow(tab.bar)){
          tab.bar[i,which(tab.bar[i,]=="NaN")]<-0
        }
      }      
      out.list<-c(out.list,list(data.frame(tab.bar,check.names=F)))
      names(out.list)<-c(names(out.list)[which(names(out.list)!="")],"gene_usage_vs_junction_frame")
      
    }
  }
  return(out.list)
}


