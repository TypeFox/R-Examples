## Julia Bischof
## 2016-02-24

plotSequencesMutationAA<-function(mutationAAtab = NULL, showChange=c("no","hydropathy","volume","chemical"), title=NULL, PDF=NULL){
  if(length(mutationAAtab)==0){
    stop("--> Output of sequences.mutation.AA() is missing")
  }
  
  if(length(showChange)!=1 || showChange=="no"){
    showChange<-NULL
  }
    
  AArange<-c("I", "V", "L", "F", "C", "M", "A", "W", 
             "G", "T", "S", "Y", "P", "H",
             "N", "D", "Q", "E", "K", "R")
  
  AAchar<-rbind(c(rep("hydrophobic", 8), rep("neutral", 6), rep("hydrophilic", 6)),
                c("large", "medium","large","very large","small","large","very small", "very large",
                  "very small","small", "very small", "very large", "small", "medium",
                  "small","small","medium","medium","large","large"),
                c("aliphatic","aliphatic","aliphatic","aromatic", "sulfur", "sulfur","aliphatic","aromatic",
                  "alipahtic","hydroxyl","hydroxyl","aromatic","aliphatic","basic",
                  "amide","acidic","amide","acidic","basic","basic"))
  rownames(AAchar)<-c("hydropathy","volume","chemical")
  colnames(AAchar)<-c("I", "V", "L", "F", "C", "M", "A", "W", 
                      "G", "T", "S", "Y", "P", "H",
                      "N", "D", "Q", "E", "K", "R")
  if(length(PDF)>0){
    pdf(paste(PDF,"_AA-mutation.pdf",sep=""),width = 9, height = 7, pointsize = 12, onefile=F)
    par(mfrow=c(1,1), mar=c(2,5,5,2),mar=c(2,5,7,2), las=1)
  }
  
  prop.tab2<-matrix("darkblue", ncol=ncol(mutationAAtab), nrow=nrow(mutationAAtab))
    for(j in 1:nrow(mutationAAtab)){
      prop.tab2[j,which(mutationAAtab[j,]<0.5)]<-"blue"
      prop.tab2[j,which(mutationAAtab[j,]<0.2)]<-"lightblue4"
      prop.tab2[j,which(mutationAAtab[j,]<0.1)]<-"deepskyblue3"
      prop.tab2[j,which(mutationAAtab[j,]<0.05)]<-"cyan"
      prop.tab2[j,which(mutationAAtab[j,]<0.01)]<-"lightblue1"
      prop.tab2[j,which(mutationAAtab[j,]==0)]<-"white"
    }
    diag(prop.tab2)<-"gray"
    prop.tab2<-t(prop.tab2)
    
    if(length(showChange)==1 && showChange=="hydropathy"){
      prop.tab3<-matrix(NA, ncol=ncol(mutationAAtab), nrow=nrow(mutationAAtab))
      for(j in 1:nrow(mutationAAtab)){
        prop.tab3[j,which(AAchar["hydropathy",]!=AAchar["hydropathy",which(colnames(AAchar)==rownames(mutationAAtab)[j])])]<-"orange"
      }
      prop.tab3<-t(prop.tab3)
    }
  
  if(length(showChange)==1 && showChange=="volume"){
    prop.tab3<-matrix(NA, ncol=ncol(mutationAAtab), nrow=nrow(mutationAAtab))
    for(j in 1:nrow(mutationAAtab)){
      prop.tab3[j,which(AAchar["volume",]!=AAchar["volume",which(colnames(AAchar)==rownames(mutationAAtab)[j])])]<-"orange"
    }
    prop.tab3<-t(prop.tab3)
  }
  
  if(length(showChange)==1 && showChange=="chemical"){
    prop.tab3<-matrix(NA, ncol=ncol(mutationAAtab), nrow=nrow(mutationAAtab))
    for(j in 1:nrow(mutationAAtab)){
      prop.tab3[j,which(AAchar["chemical",]!=AAchar["chemical",which(colnames(AAchar)==rownames(mutationAAtab)[j])])]<-"orange"
    }
    prop.tab3<-t(prop.tab3)
  }

    
    plot(1,1,col="white",xlab=" ", ylab="from", axes=F, xlim=c(1,(ncol(prop.tab2)+7)), ylim=c(1,(ncol(prop.tab2)+1)))
    for(k in 1:ncol(prop.tab2)){
      for(j in 1:nrow(prop.tab2)){
        rect(k,j,k+1,j+1, angle=0,col=prop.tab2[k,j])
      }
    }
  if(length(showChange)>0){
    for(k in 1:ncol(prop.tab3)){
      for(j in 1:nrow(prop.tab3)){
        if(mutationAAtab[j,k]>0 && !is.na(prop.tab3[k,j])){
          draw.circle(k+0.5,j+0.5,radius=0.1,col=prop.tab3[k,j])
        }
      }
    }
  }

    axis(3,at=seq(1.5,(ncol(prop.tab2)+0.5),1), AArange, tick = F, cex.axis=0.8)
    axis(2, at=seq(1.5,(ncol(prop.tab2)+0.5),1), AArange, tick = F, cex.axis=0.8)
    mtext(text = "to", side = 3, line = 3, at = ncol(prop.tab2)/2)
    mtext(text = title, side = 3, line = 5, at = ncol(prop.tab2)/2)
  if(length(showChange)>0){
    legend("topright", col=c("gray","white", "lightblue1","cyan","deepskyblue3", "lightblue4","blue","darkblue","orange"), c("same amino acid", "0%", "<1%", "1-5%", "5-10%", "10-20%", "20-50%", "50-100%", paste(showChange,"change",sep=" ")), 
           pch=c(rep(15,8),16), bg="gray90", y.intersp=1.3, cex=0.8)   
  }else{
    legend("topright", col=c("gray","white", "lightblue1","cyan","deepskyblue3", "lightblue4","blue","darkblue"), c("same amino acid", "0%", "<1%", "1-5%", "5-10%", "10-20%", "20-50%", "50-100%"), 
           pch=c(rep(15,8)), bg="gray90", y.intersp=1.3, cex=0.8)  
  }
  
  if(length(PDF)>0){
    dev.off()
  }
}
  