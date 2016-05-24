## Julia Bischof
## 10-09-2015

plotAADistribution<-function(aaDistribution.tab=NULL,plotSeqN=FALSE,colors=NULL,PDF=NULL,...){
    if(length(aaDistribution.tab)==0){
      stop("--> Amino acid distribution tab is missing or empty")
    }
    color<-c("lightblue2","red","lightgreen","lightgray","yellow","orange","purple4","cyan","darkblue","darkred",
             "darkgreen","gray33","blueviolet","aquamarine","bisque4","violetred","deepskyblue4","black","olivedrab","olivedrab3",
             "orchid1")
    if(length(PDF)>0){
      pdf(file = paste(PDF,"_Amino-acid-distribution.pdf",sep=""),
          width = ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1))*5.3,
          height = floor(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1))*6.6,
          pointsize = if(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)*3.5<32){32}else{sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)*3.5},onefile = F)
    }
    par(mfrow=if(ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1))*floor(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1))>length(aaDistribution.tab$Amino_acid_distribution)+1){c(ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)),floor(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)))}else{c(ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)),ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)))},
        mar=c(4,4,3,1),oma=c(2,2,7,3))
    for(i in 1:length(aaDistribution.tab$Amino_acid_distribution)){
      barplot(aaDistribution.tab$Amino_acid_distribution[[i]],col=color,axes=F,ylab="Percentage", cex.names = 0.9, cex.main=1.2,
              names=seq(1,ncol(aaDistribution.tab$Amino_acid_distribution[[i]]),1),xlab="Position",
              main=names(aaDistribution.tab$Amino_acid_distribution)[i])
      axis(2,at=seq(0,1,0.25),seq(0,100,25))
    }
    plot(x=seq(0.5,6.5,1),y=rep(5,7),pch=18,cex=2,col=color[1:7],axes=F,xlab="",ylab="",xlim=c(0.2,7.2),ylim=c(0,6))
    points(x=seq(0.5,6.5,1),y=rep(3,7),pch=18,cex=2,col=color[8:14])
    points(x=seq(0.5,6.5,1),y=rep(1,7),pch=18,cex=2,col=color[15:21])
    text(x=seq(1,7,1),y=5,rownames(aaDistribution.tab$Amino_acid_distribution[[1]])[1:7],cex=1)
    text(x=seq(1,7,1),y=3,rownames(aaDistribution.tab$Amino_acid_distribution[[1]])[8:14],cex=1)
    text(x=seq(1,7,1),y=1,rownames(aaDistribution.tab$Amino_acid_distribution[[1]])[15:21],cex=1)
    
    title("Amino acid distribution",outer=T, cex.main=2)
    if(length(PDF)>0){
      dev.off()
    }
  if(plotSeqN==T){    
    if(nrow(aaDistribution.tab$Number_of_sequences_per_length)==0){
      stop("--> 'Number of sequences' tab is missing or empty")
    }
    if(length(PDF)>0){
      pdf(file = paste(PDF,"_Number-of-sequences.pdf",sep=""),
          width = 7, 
          height = 9, 
          pointsize = 12) 
    }
    par(mfrow=c(1,1),mar=c(3,5,4,1))
    if(length(colors)==0){
      color<-rainbow(n=nrow(aaDistribution.tab$Number_of_sequences_per_length),start=0.1, end=0.9)
    }else{
      color<-colors
    }
    
    miny<-min(aaDistribution.tab$Number_of_sequences_per_length,na.rm=T)-min(aaDistribution.tab$Number_of_sequences_per_length,na.rm=T)/10
    maxy<-max(aaDistribution.tab$Number_of_sequences_per_length,na.rm=T)+max(aaDistribution.tab$Number_of_sequences_per_length,na.rm=T)/10
    plot(1,1,col="white",ylab="Number of sequences",xlab=" ",main="Number of sequences per length",xlim=c(0,1.7),
         ylim=c(miny,maxy),axes=F)
    axis(2,at=round(miny+(maxy-miny)*seq(0,1,0.25),0),round(miny+(maxy-miny)*seq(0,1,0.25),0))
    
    newtab<-data.frame(aaDistribution.tab$Number_of_sequences_per_length[order(aaDistribution.tab$Number_of_sequences_per_length[,1]),],
                       row.names=as.character(row.names(aaDistribution.tab$Number_of_sequences_per_length))[order(aaDistribution.tab$Number_of_sequences_per_length[,1])])
    for(i in 1:nrow(newtab)){
      lines(x = c(0,1), y = c(newtab[i,1],newtab[i,1]),col=color[i])
      text(x = 1.3, y = newtab[i,1],col=color[i], 
           labels = rownames(newtab)[i], cex=0.8)
    }
    if(length(PDF)>0){
      dev.off()
    }
  }
}


