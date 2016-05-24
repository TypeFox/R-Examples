## Julia Bischof
## 10-09-2015

plotClonesCDR3Length<-function(CDR3Length=NULL,functionality=NULL, junctionFr=NULL,
                             color=c("orange","darkblue","gray"),abundance=c("relative","absolute"),title=NULL, PDF=NULL,...){
  
  if(length(CDR3Length)==0){
    stop("--> CDR3 length vector is missing")
  }
  
  if(length(abundance)!=1 || !(abundance %in% c("relative","absolute"))){
    abundance<-"relative"
  }
  
  if(length(functionality)>0){
    functionality.new<-unlist(apply(data.frame(functionality),1,function(x){strsplit(x,split=" |,|[.]|;|[|]|_|-")[[1]]}))
    functionality.new<-functionality.new[which(nchar(functionality.new)>0)]
  }
  
  if(length(junctionFr)>0){
    junctionFr.new<-unlist(apply(data.frame(junctionFr),1,function(x){strsplit(x,split=" |,|[.]|;|[|]|_")[[1]]}))
    junctionFr.new<-junctionFr.new[which(nchar(junctionFr.new)>0)]
  }
  

    minCDR3<-min(as.numeric(CDR3Length),na.rm=T)
    maxCDR3<-max(as.numeric(CDR3Length),na.rm=T)
    # only CDR3 length
      tab.bar<-t(data.frame(apply(data.frame(seq(minCDR3,maxCDR3,1)),1,function(x){length(which(CDR3Length==x))})))
      colnames(tab.bar)<-as.character(seq(minCDR3,maxCDR3,1))
      if(abundance=="relative"){
        tab.bar<-tab.bar/length(CDR3Length)*100
      }
    if(length(PDF)>0){
      pdf(file = paste(PDF,"_CDR3-length.pdf",sep=""),width = if(ncol(tab.bar)*0.7<7){7}else{ncol(tab.bar)*0.7},height = 7,pointsize = if(ncol(tab.bar)*0.5<16){16}else{ncol(tab.bar)*0.5})
    }
      barplot(tab.bar,col=color,xlab="CDR3 length",ylab=if(abundance=="relative"){"Percentage"}else{"Quantity"},main=title, cex.names = 0.8,
              ylim=c(0,max(tab.bar)+(max(tab.bar)/5)))
  
  
  if(length(PDF)>0){
      dev.off()
    }
    # CDR3 length vs. Functionality
    if(length(CDR3Length)==length(functionality)){
      CDR3Length.new<-vector()
      for(i in 1:length(CDR3Length)){
        CDR3Length.new<-c(CDR3Length.new,rep(CDR3Length[i],length(grep(",",functionality[i]))+1))
      }
      tab.bar<-vector()
      for(i in c("^productive","^unproductive","unknown|No results")){
        tab.bar<-rbind(tab.bar,apply(data.frame(seq(minCDR3,maxCDR3,1)),1,function(x){length(intersect(grep(i,functionality.new),which(CDR3Length.new==x)))}))
      }
      colnames(tab.bar)<-as.character(seq(minCDR3,maxCDR3,1))
      rownames(tab.bar)<-c("productive","unproductive","unknown")
      if(abundance=="relative"){
        tab.bar<-t(data.frame(apply(tab.bar,1,function(x){x/colSums(tab.bar,na.rm=T)*100})))
        for(i in 1:nrow(tab.bar)){
          tab.bar[i,which(tab.bar[i,]=="NaN")]<-0
        }
      }
      if(length(PDF)>0){
        pdf(file = paste(PDF,"_CDR3-length_vs_functionality.pdf",sep=""),width = if(ncol(tab.bar)*0.7<7){7}else{ncol(tab.bar)*0.7},height = 7,pointsize = if(ncol(tab.bar)*0.5<16){16}else{ncol(tab.bar)*0.5})
      }
      barplot(tab.bar,col=color,xlab="CDR3 length",ylim=if(abundance=="absolute"){c(0,max(tab.bar)+(max(tab.bar)/5))}else{c(0,160)},main=title,
              ylab=if(abundance=="relative"){""}else{"Quantity"},axes=if(abundance=="relative"){F}else{T}, cex.names = 0.8)
      if(abundance=="relative"){
        axis(2,at=seq(0,100,25),seq(0,100,25))
        mtext(side = 2,at = 50,line = 2.5,text = "Percentage")
      }
      legend("topright",col=color,c("productive","unproductive","unknown"),title="Functionality",pch=15,y.intersp=1.1)
      if(length(PDF)>0){
        dev.off()
      }
    }else if(length(functionality)>0 && length(CDR3Length)!=length(functionality)){
      stop("--> Lengths of CDR3 Length vector and functionality vector are different")
    }
    
    # CDR3 length vs. JUNCTION_frame
    if(length(CDR3Length)==length(junctionFr)){
      CDR3Length.new<-vector()
      for(i in 1:length(CDR3Length)){
        CDR3Length.new<-c(CDR3Length.new,rep(CDR3Length[i],length(grep(",",junctionFr[i]))+1))
      }
      tab.bar<-vector()
      for(i in c("in-frame","out-of-frame","null")){
        tab.bar<-rbind(tab.bar,apply(data.frame(seq(minCDR3,maxCDR3,1)),1,function(x){length(intersect(grep(i,junctionFr.new),which(CDR3Length.new==x)))}))
      }
      colnames(tab.bar)<-as.character(seq(minCDR3,maxCDR3,1))
      rownames(tab.bar)<-c("in-frame","out-of-frame","unknown")
      if(abundance=="relative"){
        tab.bar<-t(data.frame(apply(tab.bar,1,function(x){x/colSums(tab.bar,na.rm=T)*100})))
        for(i in 1:nrow(tab.bar)){
          tab.bar[i,which(tab.bar[i,]=="NaN")]<-0
        }
      }
        if(length(PDF)>0){
          pdf(file = paste(PDF,"_CDR3-length_vs_junction_frame.pdf",sep=""),width = if(ncol(tab.bar)*0.7<7){7}else{ncol(tab.bar)*0.7},height = 7,pointsize = if(ncol(tab.bar)*0.5<16){16}else{ncol(tab.bar)*0.5})
        }
      barplot(tab.bar,col=color,xlab="CDR3 length",ylim=if(abundance=="absolute"){c(0,max(tab.bar)+(max(tab.bar)/5))}else{c(0,160)},main=title,
              ylab=if(abundance=="relative"){""}else{"Quantity"},axes=if(abundance=="relative"){F}else{T}, cex.names = 0.8)
      if(abundance=="relative"){
        axis(2,at=seq(0,100,25),seq(0,100,25))
        mtext(side = 2,at = 50,line = 2.5,text = "Percentage")
      }
      legend("topright",col=color,c("in-frame","out-of-frame","unknown"),title="Junction frame",pch=15,y.intersp=1.1)
      if(length(PDF)>0){
        dev.off()
      }
    }else if(length(junctionFr)>0 && length(CDR3Length)!=length(junctionFr)){
      stop("--> Lengths of CDR3 Length vector and junction frame vector are different")
    }
}


