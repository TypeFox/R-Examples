plotCluster <- function(DRdata,doseData,ORCMEoutput,clusterID,zeroMean=FALSE,xlabel,ylabel,main=""){
	plotData <- DRdata[ORCMEoutput[, clusterID], ]
	if (zeroMean) {
		plotData2 <- plotData - rowMeans(plotData)
	}
	else {
		plotData2 <- plotData
	}
    ngene <- dim(plotData2)[1]
      
	plot(as.numeric(sort(unique(doseData))), plotData2[1, ], ylim = range(plotData2), 
			xlab = xlabel, ylab = ylabel, cex.axis = 1.5, cex = 1.5, 
			cex.lab = 1.5, type = "l", xaxt = "n")
	axis(1, at = c(1:length(unique(doseData))), labels = sort(unique(doseData)), 
			cex.axis = 1.5, cex = 1.5, cex.lab = 1.5)
	apply(plotData2, 1, function(x) lines(sort(unique(doseData)), x, 
						cex.axis = 1.5, cex = 1.5, cex.lab = 1.5))
}


