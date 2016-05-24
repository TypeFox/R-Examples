
optimalGT <- function(delta.t.mat.allmk.list,gt0, gt.thres=0,optGTplot=FALSE){
  #input  delta.t.mat.allmk.list
  #a list with length=nMK, each element is a matrix nSample by nSigGene ,
  # each element matrix is deltaT for each sigGene(row) when  sample i (column) is flipped
  #gt0 is the original genotype
  #when mean deltaT of all top sig traits at certain marker is > gt.thres, the old gt is flipped (0-->1 or 1-->0)
  
  
           optGT <- gt0
           wls.mat <- matrix(0, nrow(gt0), ncol(gt0)) #nMk x nSample
           colnames(wls.mat) <- colnames(gt0)
           for(i.mk in 1:length(delta.t.mat.allmk.list)){
              delta.t.mat   <-      delta.t.mat.allmk.list[[i.mk]]  #a matrix nSample by nSigGene
              aveDeltaT     <- apply(delta.t.mat,2, function(x) median(x,na.rm=TRUE))
              temp          <- abs(gt0[i.mk,which(aveDeltaT>gt.thres)] -1)           #check the column of delta.t.mat is same order as gt data
              #cat(i.mk,"\n", which(aveDeltaT>gt.thres),"\n", file="temp.txt",append=T)
              optGT[i.mk,which(aveDeltaT>gt.thres)] <- temp
              wls.mat [i.mk,which(aveDeltaT>gt.thres)] <- 1
           }
           if(optGTplot) {
            #library('gplots')
            kleuren=c("green","red")
            mymain="" ;plotName="WLS_mat_01plot"
            png(paste(plotName, ".png",sep=""))
            heatmap.2(wls.mat,dendrogram="none", key=F, Rowv =F,scale="none",main=mymain,trace="none",colsep=seq(1,ncol(gt0),by=1),density.info="none",col=kleuren,breaks<-seq(-0.5,1.2,length=3)  )
            dev.off()
           }
          return(optGT)
}
