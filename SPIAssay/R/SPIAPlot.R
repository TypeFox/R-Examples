############################################################################
#
# SPIAPlot
#  Input: 1. SPIAanalysis: the result of SPIATest function
#  Output: a plot of the cell line pairs versus the SPIA distance with 
#          information about the probabilistic test
#
#############################################################################

SPIAPlot <- function(SPIAanalysis){  
  
  #set parameters
  par(mfrow=c(1,1),cex=0.8,font=1,col.main="skyblue4",mai=c(1, 0.8, 0.3, 0.1))
  titleplot<-paste("Pair-wise comparison of ", SPIAanalysis$input.param$N_samples," samples on ", SPIAanalysis$input.param$N_SNPs, " SNPs",sep="")  
  plot(1:10,xlim=c(0,dim(SPIAanalysis$SPIAresult)[1]),ylim=c(0,1),ylab="Distance D (% of discordant genotype calls)",xlab="Index of Cell Line pairs (all combinations)",main=titleplot,col="white" ,cex=0.1);  

  #check if the SPIAanalysis contain the probabilistic test
  if (SPIAanalysis$input.param$testDone)
  #compute the color depending on the result of the statistical test
  {      
    for(i in c(1:dim(SPIAanalysis$SPIAresult)[1])){
      pointConf <- switch(SPIAanalysis$SPIAresult[i,4], Uncertain = c("blue",19), Similar = c("springgreen3",19), Different = c("red",20), c("black",19))    
      points(i,SPIAanalysis$SPIAresult[i,3],col=pointConf[1],pch=as.integer(pointConf[2]),cex=0.7);     
    }
    #plot the legend that describe the color code
    legendText <- c("SPIA TEST: different","SPIA TEST: uncertain","SPIA TEST: match", paste("< ",SPIAanalysis$parameters$PercValidCall,"% of available calls",sep=""))
    legend(0,1,legendText, col = c("red","blue","springgreen3","black"), pch=c(20,19,19,19),cex=0.6)                            
  } else {
    #use the same color for each point
    for(i in c(1:dim(SPIAanalysis$SPIAresult)[1])){      
      points(i,SPIAanalysis$SPIAresult[i,3],col="black",pch=20,cex=0.7);     
    }
  }
  
  
  

}