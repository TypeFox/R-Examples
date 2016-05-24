HeatMapFactor <-
function(pts,class,factor,SaveFileName="",Title=""){
  # a heatmap figure for visualizing how a SNP or a SNP-combination influencing
  # the phenotype.
  #
  # input
  #    pts: Row -> Sample, Column -> SNP
  #         1 -> AA
  #         2 -> Aa
  #         3 -> aa
  #    class: Row -> 1, Column -> class label
  #           1 -> case
  #           2 -> control
  #    factor: the considered SNP or SNP-combination, for example, factor <- 5
  #            or factor <- c(2,5)
  #    SaveFileName: basic file name for saving figure.
  #    Title: Tile of the figure
  #
  # Junliang Shang
  # 4.24/2014
  
  #library(ggplot2)
  #library(reshape2)
  
  ## data
  DataNum <- length(factor)
  
  SNPdata <- rep(0,3^DataNum)
  Totaldata <- rep(0,3^DataNum)
  
  for (i in 1:nrow(pts)){
    info <- 0;
    count <- 0;
    
    for (j in 1:DataNum){
      if (pts[i,factor[j]]!=0){
        count <- (pts[i,factor[j]]-1)*(3^(DataNum-j))+count;
      }else
      {
        info <- 1
        break
      }
    }
    
    if (info==0){
      Totaldata[count+1] <- Totaldata[count+1]+1
      
      if (class[i]==1){
        SNPdata[count+1] <- SNPdata[count+1]+1
      }
    }
  }
  
  SNPdata <- SNPdata/Totaldata
  
  SNPdata[is.nan(SNPdata)] <- 0

  # translation
  
  if(DataNum==1){
    SNPdata <- t(matrix(SNPdata,1,3))
    SNPdata <- as.data.frame(SNPdata)
    SNPdata <- cbind(c("AA","Aa","aa"), SNPdata)
    colnames(SNPdata) <- c("Genotype","")
  }
   
  if(DataNum==2){
    SNPdata <- t(matrix(SNPdata,3,3))
    SNPdata <- as.data.frame(SNPdata)
    SNPdata <- cbind(c("AA","Aa","aa"), SNPdata)
    colnames(SNPdata) <- c("Genotype","BB","Bb","bb")
  }
  
  if(DataNum==3){
    SNPdata <- t(matrix(SNPdata,9,3))
    SNPdata <- as.data.frame(SNPdata)
    SNPdata <- cbind(c("AA","Aa","aa"), SNPdata)
    colnames(SNPdata) <- c("Genotype","BBCC","BBCc","BBcc","BbCC","BbCc","Bbcc","bbCC","bbCc","bbcc")
  }
   
  if(DataNum==4){
    SNPdata <- t(matrix(SNPdata,9,9))
    SNPdata <- as.data.frame(SNPdata)
    SNPdata <- cbind(c("AABB","AABb","AAbb","AaBB","AaBb","Aabb","aaBB","aaBb","aabb"), SNPdata)
    colnames(SNPdata) <- c("Genotype","CCDD","CCDd","CCdd","CcDD","CcDd","Ccdd","ccDD","ccDd","ccdd")
  }
  
  if(DataNum==5){
    SNPdata <- t(matrix(SNPdata,27,9))
    SNPdata <- as.data.frame(SNPdata)
    SNPdata <- cbind(c("AABB","AABb","AAbb","AaBB","AaBb","Aabb","aaBB","aaBb","aabb"), SNPdata)
    colnames(SNPdata) <- c("Genotype","CCDDEE","CCDDEe","CCDDee","CCDdEE","CCDdEe","CCDdee","CCddEE","CCddEe","CCddee",
                           "CcDDEE","CcDDEe","CcDDee","CcDdEE","CcDdEe","CcDdee","CcddEE","CcddEe","Ccddee",
                           "ccDDEE","ccDDEe","ccDDee","ccDdEE","ccDdEe","ccDdee","ccddEE","ccddEe","ccddee")
  }
  
  if(!(DataNum %in% c(1,2,3,4,5))){
    stop(" MaxOrder Error!\n")
  }
  
  results <- SNPdata
  SNPdata <- melt(SNPdata,id.vars="Genotype")
  
  # print
  dev.new()
  set.seed(8)
  
  variable <- Genotype <- value <- NULL
  
  p <- ggplot(SNPdata, aes(variable, Genotype))+xlab(" ")+ylab(" ")+labs(title=Title)
  p <- p + geom_tile(aes(fill = value), colour = "white")
  p <- p+scale_fill_gradient(low='green', high='red')
  p <- p+geom_text(aes(label=round(value,3)), angle=45)
  
  print(p)
  
  # save and return
  # savePlot(filename=paste(SaveFileName,"_HF",sep=""),type="tiff")
  list(HeatMapFactors=results)
}
