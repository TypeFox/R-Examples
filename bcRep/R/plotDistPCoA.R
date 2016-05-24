# Julia Bischof
# 12-08-2015

plotDistPCoA<-function(pcoa.tab=NULL, groups=NULL, names=NULL, axes=NULL, plotCorrection=FALSE, title=NULL, plotLegend=FALSE, PDF=NULL){
  if(length(pcoa.tab)==0){
    stop("--> PCoA matrix is missing")
  }
  if(length(axes)==0){
    axes<-c(1,2)
  }
  if(length(axes)!=2){
    stop("--> Insert two axes: f.e. c(1,2)")
  }
  
  minaxes<-Inf
  for(i in 1:length(pcoa.tab)){
    if(pcoa.tab[[i]][1]!="no results available" && ncol(pcoa.tab[[i]]$vectors)<minaxes){
      minaxes<-ncol(pcoa.tab[[i]]$vectors)
    }
  }
  if(minaxes<axes[1] || minaxes<axes[2]){
    stop(paste("--> Not enough axes available (minimum number of axes = ",minaxes,")",sep=""))
  }
  if(length(groups)==0){
    stop("--> 'groups' has to be a data.frame containing sequences (1. column) and groups (2.column)")  }
  
  if(length(title)==0){
    title2<-"Principal coordinate analysis"
  }else{
    title2<-title
  }
  
  
  if(nrow(groups)>0 && ncol(groups)==2){
    temp<-rainbow(n=length(unique(groups[,2])),start=0.2,end=0.9)
    colors<-vector("character", length = nrow(groups))
    for(c in 1:length(unique(groups[,2]))){
      colors[which(groups[,2]==unique(groups[,2])[c])]<-temp[c]
    }
  }else{
    stop("--> 'groups' has to be a data.frame containing sequences (1. column) and groups (2.column)")
  }
  
  if(length(PDF)>0){
    pdf(paste(PDF,"_sequences_PCoA.pdf",sep=""),width = if(length(pcoa.tab)==1){10}else{ceiling(sqrt(length(pcoa.tab)+1))*4}, 
        height = if(length(pcoa.tab)==1){5}else{round(sqrt(length(pcoa.tab)+1))*4}, 
        pointsize =if(length(pcoa.tab)==1){12}else{ceiling(sqrt(length(pcoa.tab)+1))*4})
  }
  if(length(pcoa.tab)>1 && plotLegend==T){
    par(mfrow=c(round(sqrt(length(pcoa.tab)+1)),ceiling(sqrt(length(pcoa.tab)+1))), mar=c(5,5,5,2), oma=c(1,1,6,1))
  }else if(length(pcoa.tab)>1 && plotLegend==F){
    par(mfrow=c(round(sqrt(length(pcoa.tab))),ceiling(sqrt(length(pcoa.tab)))), mar=c(5,5,5,2), oma=c(1,1,6,1))
  }else if(length(pcoa.tab)==1 && plotLegend==T){
    par(mfrow=c(1,2),mar=c(5,5,0,1),oma=c(1,1,6,1))
  }else if(length(pcoa.tab)==1 && plotLegend==F){
    par(mfrow=c(1,1),mar=c(5,5,0,1),oma=c(1,1,6,1))
  }
  for(i in 1:length(pcoa.tab)){
    color.vec<-vector()
    if(as.character(pcoa.tab[[i]][1])=="no results available"){
      plot(1,1,type="n", xlab="",ylab="",axes=F, main=names(pcoa.tab)[i])
      text(1,1,"no results available")
    }else{
        temp<-cbind(groups, colors)
        temp<-temp[which(temp[,1] %in% rownames(pcoa.tab[[i]]$vectors)),]
        color.vec<-temp[,3]
        if(plotCorrection==TRUE){
          plot(pcoa.tab[[i]]$vectors.cor[,axes[1]],pcoa.tab[[i]]$vectors.cor[,axes[2]], 
               xlab=paste("axis ", axes[1], " (",round(pcoa.tab[[i]]$values$Rel_corr_eig[axes[1]]/sum(pcoa.tab[[i]]$values$Rel_corr_eig,na.rm=T)*100,1),"%)",sep=""),
               ylab=paste("axis ", axes[2], " (",round(pcoa.tab[[i]]$values$Rel_corr_eig[axes[2]]/sum(pcoa.tab[[i]]$values$Rel_corr_eig,na.rm=T)*100,1),"%)",sep=""), 
               main=names(pcoa.tab)[i], col=color.vec, pch=19)
        }else{
          plot(pcoa.tab[[i]]$vectors[,axes[1]],pcoa.tab[[i]]$vectors[,axes[2]], 
               xlab=paste("axis ", axes[1], " (",round(pcoa.tab[[i]]$values$Rel_corr_eig[axes[1]]/sum(pcoa.tab[[i]]$values$Rel_corr_eig,na.rm=T)*100,1),"%)",sep=""),
               ylab=paste("axis ", axes[2], " (",round(pcoa.tab[[i]]$values$Rel_corr_eig[axes[2]]/sum(pcoa.tab[[i]]$values$Rel_corr_eig,na.rm=T)*100,1),"%)",sep=""), 
               main=names(pcoa.tab)[i], col=color.vec, pch=19)
        }      
    }
  }
  if(plotLegend==T){
    plot(1,1,xlab="", ylab="", axes=F, type="n")
    legend("center", col=unique(colors), unique(groups[,2]), pch=19, y.intersp = 1.1)
  }
  title(title2,outer=T, cex.main=1.5)
  
  if(length(PDF)>0){
    dev.off()
  }
  
}
