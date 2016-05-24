#' plotcsSAM
#' 
#' Plots the # of genes called significnat at a given false disocvery rate for
#' the SAM (heterogenous tissue) comparison, and for each of the contrasted
#' cell-types using csSAM
#' 
#' 
#' @param csSAMdata List object output of the fdrCsSAM function
#' @param SAMdata List object output of the fdrSAM function
#' @param alternative Type of test conducted. Will appear in plot title.
#' @param cellID Label for each cell-type
#' @param numcell Number of different cell-types being considered.
#' @param fileName Name of output pdf file.
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @cite Shen-Orr2010
plotCsSAM <-
function(csSAMdata, SAMdata,alternative,cellID,numcell,fileName){
  
  cellnames = cellID
  labs=paste("cell type ",as.character(1:numcell),sep="")
  
  pdf(file=fileName, height=11,width=8)
   
  par(mfrow=c(3,3))
  plot(SAMdata$ncall.sam,SAMdata$fdr.sam,xlab="## called",ylab="FDR", type="l",log="x",ylim=c(0,1))
  title(paste("SAM",alternative))
  for(i in 1:numcell) {
    plot(csSAMdata$ncall.g[i,],csSAMdata$fdr.g[i,],xlab="# called",ylab="FDR", type="l",log="x",ylim=c(0,1))
    title(paste(cellnames[i],alternative))
  }
  dev.off()
}

