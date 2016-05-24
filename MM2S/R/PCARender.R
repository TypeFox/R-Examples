# Deena M.A. Gendoo
# October 29, 2014
# Code to conduct PCA projections of MM2S predicted (versus original designations) of MB subtypes

# DISCLAIMER:
# MM2S package (and its code components) is provided "AS-IS" and without any warranty of any kind. 
# In no event shall the University Health Network (UHN) or the authors be liable for any consequential damage of any kind, 
# or any damages resulting from the use of this MM2S.

#################################################################################
#################################################################################

PCARender<-function(GSVAmatrixTesting,GSVAmatrixTraining)
{
  message("Three PDFs have been generated, please consult your working directory to find them.")
  
  Training<-GSVAmatrixTraining
  col.list<-c("goldenrod1","#abdda4","#bababa","#d7191c","#2b83ba")[as.factor(MB_SampleInfo$subtype)]
  pch.list<-c(rep(16,nrow(Training)))
  cex.list<-c(rep(0.75,nrow(Training)))
  
  pcaTraining<-prcomp(Training,scale=T)
  SpikeTestOntoMB<-predict(pcaTraining,newdata=GSVAmatrixTesting)
  FullMatrix<-rbind(pcaTraining$x[,1:3],SpikeTestOntoMB[,1:3])

  pdf("PC1-PC2.pdf")
  plot(pcaTraining$x[,1:2], type="p",pch=pch.list, cex=cex.list,col=col.list,xlim=c(min(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2]), max(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2])),ylim=c(min(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2]), max(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2])))
  par(new=T)
  points(SpikeTestOntoMB[,1:2], type="p",pch=18, cex=1,col="purple",xlim=c(min(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2]), max(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2])),ylim=c(min(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2]), max(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2])))
  text(SpikeTestOntoMB[,1:2],rownames(SpikeTestOntoMB),cex=0.4)
  legend("bottom",fill=c("goldenrod1","#abdda4","#bababa","#d7191c","#2b83ba"),legend=c("Group3","Group4","Normal","SHH","WNT"),horiz=T,cex=0.5)
  dev.off()
  
  pdf("PC2-PC3.pdf")
  plot(pcaTraining$x[,2:3], type="p",pch=pch.list, cex=cex.list,col=col.list,xlim=c(min(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3]), max(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3])),ylim=c(min(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3]), max(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3])))
  par(new=T)
  points(SpikeTestOntoMB[,2:3], type="p",pch=18, cex=1,col="purple",xlim=c(min(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3]), max(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3])),ylim=c(min(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3]), max(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3])))  
  text(SpikeTestOntoMB[,2:3],rownames(SpikeTestOntoMB),cex=0.4)
  legend("bottom",fill=c("goldenrod1","#abdda4","#bababa","#d7191c","#2b83ba"),legend=c("Group3","Group4","Normal","SHH","WNT"),horiz=T,cex=0.5)
  dev.off()

  col.list2<-c(rep("purple",nrow(SpikeTestOntoMB)))
  col.full<-c(col.list,col.list2)
  pch.list2<-c(rep(17,nrow(SpikeTestOntoMB)))
  pch.full<-c(pch.list,pch.list2)
  cex.full<-c(rep(0.55,nrow(FullMatrix)))
  pdf("LatticeMatrix.pdf")
  print(splom(as.data.frame(FullMatrix),type="p",pch=pch.full,col=col.full,xlim=c(min(FullMatrix), max(FullMatrix)), ylim=c(min(FullMatrix), max(FullMatrix)),cex=cex.full))
  dev.off()
  
  #First Plot PC1-PC2
  plot(pcaTraining$x[,1:2], type="p",pch=pch.list, cex=cex.list,col=col.list,xlim=c(min(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2]), max(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2])),ylim=c(min(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2]), max(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2])))
  par(new=T)
  points(SpikeTestOntoMB[,1:2], type="p",pch=18, cex=1,col="purple",xlim=c(min(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2]), max(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2])),ylim=c(min(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2]), max(pcaTraining$x[,1:2], SpikeTestOntoMB[,1:2])))
  text(SpikeTestOntoMB[,1:2],rownames(SpikeTestOntoMB),cex=0.4)
  legend("bottom",fill=c("goldenrod1","#abdda4","#bababa","#d7191c","#2b83ba"),legend=c("Group3","Group4","Normal","SHH","WNT"),horiz=T,cex=0.5)
  #Second Plot
  plot(pcaTraining$x[,2:3], type="p",pch=pch.list, cex=cex.list,col=col.list,xlim=c(min(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3]), max(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3])),ylim=c(min(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3]), max(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3])))
  par(new=T)
  points(SpikeTestOntoMB[,2:3], type="p",pch=18, cex=1,col="purple",xlim=c(min(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3]), max(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3])),ylim=c(min(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3]), max(pcaTraining$x[,2:3], SpikeTestOntoMB[,2:3])))  
  text(SpikeTestOntoMB[,2:3],rownames(SpikeTestOntoMB),cex=0.4)
  legend("bottom",fill=c("goldenrod1","#abdda4","#bababa","#d7191c","#2b83ba"),legend=c("Group3","Group4","Normal","SHH","WNT"),horiz=T,cex=0.5)
  #Third Plot
  print(splom(as.data.frame(FullMatrix),type="p",pch=pch.full,col=col.full,xlim=c(min(FullMatrix), max(FullMatrix)), ylim=c(min(FullMatrix), max(FullMatrix)),cex=cex.full))
  
  
}



