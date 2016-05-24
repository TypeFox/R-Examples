## Julia Bischof
## 10-15-2015

plotCompareAADistribution<-function(comp.tab=NULL, plotSeqN=FALSE, colors=NULL, title=NULL, PDF=NULL){
  if(length(comp.tab)==0){
    stop("--> Comparison table is missing")
  }
  if(length(title)>0){
    plottitle<-title
  }else{
    plottitle<-NULL
  }
  
  uniAA<-c("F","L","I","M","V","S","P","T","A","Y","H","C","W","N","D","G","Q","R","K","E","*")
    
  uniseqlength<-vector()
  for(i in 1:length(comp.tab$Amino_acid_distribution)){
    uniseqlength<-c(uniseqlength,names(comp.tab$Amino_acid_distribution[[i]]))
  }
  uniseqlength<-unique(uniseqlength)
  orderCDR3<-as.numeric(unlist(apply(data.frame(uniseqlength),1,function(x){strsplit(x, split=" = ")[[1]][2]})))
  uniseqlength<-uniseqlength[order(orderCDR3)]
  
  color<-c("lightblue2","red","lightgreen","lightgray","yellow","orange","purple4","cyan","darkblue","darkred",
           "darkgreen","gray33","blueviolet","aquamarine","bisque4","violetred","deepskyblue4","black","olivedrab","olivedrab3",
           "orchid1")
  if(length(PDF)>0){
    pdf(file = paste(PDF,"_Comparison_Amino-acid-distribution.pdf",sep=""),
        width = if(length(comp.tab$Amino_acid_distribution)*2<8){10}else{length(comp.tab$Amino_acid_distribution)*2+3},
        height = if(length(uniseqlength)<7){9}else{length(uniseqlength)+5},
        pointsize = if(length(uniseqlength)<10){16}else{length(uniseqlength)*0.4})
  }
  
  par(mfrow=c(length(uniseqlength)+2,length(comp.tab$Amino_acid_distribution)+1),
      mar=c(0.5,0.5,0.5,0.5),oma=c(2,2,6,2))
  plot(0,0,col="white",xlab="",ylab="",axes=F)
  for(i in 1:length(comp.tab$Amino_acid_distribution)){
    plot(1,1,col="white",xlab="",ylab="",axes=F)
    text(1,1,names(comp.tab$Amino_acid_distribution)[i],cex = 1.3)
  }
    
  for(i in 1:length(uniseqlength)){    
    plot(1,1,col="white",xlab="",ylab="",axes=F)
    text(1,1,paste("Length = ",strsplit(uniseqlength[i],split=" ")[[1]][length(strsplit(uniseqlength[i],split=" ")[[1]])], " AA",sep=""),cex=1.3)
    
    for(j in 1:length(comp.tab$Amino_acid_distribution)){
      if(uniseqlength[i] %in% names(comp.tab$Amino_acid_distribution[[j]])){
        xx<-barplot(as.matrix(as.data.frame(comp.tab$Amino_acid_distribution[[j]][which(names(comp.tab$Amino_acid_distribution[[j]])==uniseqlength[i])])),
                ylim=c(0,1),axes=F,col=color,names=rep("",as.numeric(strsplit(uniseqlength[i],split=" ")[[1]][length(strsplit(uniseqlength[i],split=" ")[[1]])])))
      }else{
        plot(1,1,col="white",xlab="",ylab="",axes=F)
      }
    }
  }
  if((length(comp.tab$Amino_acid_distribution)+1)<4){
    plot(1,1,col="white",xlab="",ylab="",axes=F)
    for(i in 2:3){
      plot(1,1,col="white",xlab="",ylab="",axes=F)
      legend("top",legend = uniAA[((5*(i-1))-4):(5*(i-1))], col=color[((5*(i-1))-4):(5*(i-1))], pch=18, horiz = TRUE,cex=0.85,pt.cex=1.8, box.col = "white")
      legend("bottom",legend = uniAA[((5*i)-4):(5*i)], col=color[((5*i)-4):(5*i)], pch=18, horiz = TRUE,cex=0.85,pt.cex=1.8, box.col = "white")
      
    }
  }else{
    for(i in 1:4){
      plot(1,1,col="white",xlab="",ylab="",axes=F)
      legend("right",legend = uniAA[((5*i)-4):(5*i)], col=color[((5*i)-4):(5*i)], pch=18, horiz = TRUE,cex=0.85,pt.cex=1.8, box.col = "white")
    }
  }

  
  title(if(length(plottitle)>0){plottitle}else{"Amino acid distribution"},outer=T, cex.main=2)
  if(length(PDF)>0){
    dev.off()
  }
  if(plotSeqN==T){    
    if(length(comp.tab$Number_of_sequences_per_length)==0){
      stop("--> 'Number of sequences' tab is missing or empty")
    }
    
    nrleg<-ceiling(length(comp.tab$Number_of_sequences_per_length)/4)
    
    if(length(PDF)>0){
      pdf(file = paste(PDF,"_Comparison_Number-of-sequences.pdf",sep=""),
          width = sqrt((length(uniseqlength)+nrleg+1))*2.3,
          height = sqrt((length(uniseqlength)+nrleg+1))*1.5, 
          pointsize = if(sqrt(length(uniseqlength)+nrleg)*2.1<12){12}else{sqrt(length(uniseqlength)+nrleg)*2.1})
    }
    if(length(colors)!=length(comp.tab$Number_of_sequences_per_length)){
      color<-rainbow(n=length(comp.tab$Number_of_sequences_per_length),start = 0.2, end = 0.8)
    }else{
      color<-colors
    }
    par(mfrow=c(ceiling(sqrt((length(uniseqlength)+nrleg+1))),ceiling(sqrt((length(uniseqlength)+nrleg+1)))),
        mar=c(1,3,1,1),oma=c(2,2,6,3))
    for(i in 1:length(uniseqlength)){
      maxy<-0
      miny<-Inf
      for(j in 1:length(comp.tab$Number_of_sequences_per_length)){
        if(length(as.numeric(comp.tab$Number_of_sequences_per_length[[j]][which(names(comp.tab$Number_of_sequences_per_length[[j]])==uniseqlength[i])]))==1 && as.numeric(comp.tab$Number_of_sequences_per_length[[j]][which(names(comp.tab$Number_of_sequences_per_length[[j]])==uniseqlength[i])])>maxy){
          maxy<-as.numeric(comp.tab$Number_of_sequences_per_length[[j]][which(names(comp.tab$Number_of_sequences_per_length[[j]])==uniseqlength[i])])
        }
        if(length(as.numeric(comp.tab$Number_of_sequences_per_length[[j]][which(names(comp.tab$Number_of_sequences_per_length[[j]])==uniseqlength[i])]))==1 && as.numeric(comp.tab$Number_of_sequences_per_length[[j]][which(names(comp.tab$Number_of_sequences_per_length[[j]])==uniseqlength[i])])<miny){
          miny<-as.numeric(comp.tab$Number_of_sequences_per_length[[j]][which(names(comp.tab$Number_of_sequences_per_length[[j]])==uniseqlength[i])])
        }
      }
      maxy<-round(maxy+maxy/10,0)
      miny<-round(miny-miny/10,0)
      
      plot(1,1,col="white",xlab=" ",main=uniseqlength[i],
           ylim=c(miny,maxy),axes=F,cex=0.8)
      if(maxy<=5){
        axis(2,at=seq(miny,maxy,1),seq(miny,maxy,1))
      }else{
        axis(2,at=c(miny,round(maxy/2,0),maxy),c(miny,round(maxy/2,0),maxy))
      }
      for(j in 1:length(comp.tab$Number_of_sequences_per_length)){
        if(uniseqlength[i] %in% names(comp.tab$Number_of_sequences_per_length[[j]])){
          abline(h=as.numeric(comp.tab$Number_of_sequences_per_length[[j]][which(names(comp.tab$Number_of_sequences_per_length[[j]])==uniseqlength[i])]),
                 col=color[j], lwd=2)
        }
      }
    }
    if(nrleg>1){
      for(j in 1:nrleg){
        plot(1,1,col="white",axes=F,xlab="",ylab="")
        legend("center",col=color[((j*4)-3):(j*4)],names(comp.tab$Number_of_sequences_per_length)[((j*4)-3):(j*4)],lty=1, lwd=3, cex=1.1, y.intersp = 1.2, box.col = "white")
      }
    }else{
      plot(1,1,col="white",axes=F,xlab="",ylab="")
      legend("center",col=color,names(comp.tab$Number_of_sequences_per_length),lty=1, lwd=3, cex=1, y.intersp = 1.1, box.col = "white")
    }
    title(if(length(plottitle)>0){plottitle}else{"Number of sequences per length"},outer=T, cex.main=2)
    if(length(PDF)>0){
      dev.off()
    }
  }
}
