plotIsomeans <- function(monoData,obsData,doseData,geneIndex){
      obsMean <- obsData[geneIndex,]
      isoMean <- monoData[geneIndex,]

      plot(doseData,obsMean,ylim=range(obsMean),xlab="Dose",ylab="Gene Expression",cex.axis=1.5,cex=1.5,cex.lab=1.5,type="p",xaxt="n",pch=19)
      axis(1, at = c(1:length(unique(doseData))),labels = c(1:length(unique(doseData))),cex.axis=1.5,cex=1.5,cex.lab=1.5)
      lines(unique(doseData),isoMean,cex.axis=1.5,cex=1.5,cex.lab=1.5,lwd=1.5)
      
}
