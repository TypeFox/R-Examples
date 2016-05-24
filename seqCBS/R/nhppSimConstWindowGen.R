nhppSimConstWindowGen <-
function(controlRates, filename, chromosomeN, nSpike=25, cptLen=c(3,5,8,12,20,30,50,75,100), nPair=2, nRepeat=10, minGain=1.5, maxGain=4, minLoss=0.01, maxLoss=0.5, pGain=0.6) {
	nCptLen = length(cptLen)
	segseqRunName = paste("segTime = zeros(", nCptLen, ",", nRepeat, ");", sep="")
	for(i in 1:nCptLen) {
		for(j in 1:nRepeat) {
			segseqRunName = c(segseqRunName, "jsStart=tic;", paste("SegSeq -i ", filename, "_Meta_", cptLen[i], "_", j, ".txt -s ",filename,cptLen[i],"_",j, ";", sep=""), paste("segTime(", i, ",", j, ")=toc(jsStart);", sep=""))
		}
	}
	segseqRunName = c(segseqRunName, paste("dlmwrite(", "'", filename, "_Time.txt'", ", segTime);", sep=""), "exit;")
	segseqRunName = as.matrix(segseqRunName, ncol=1)
	write.table(segseqRunName, file=paste(filename, "Run.m", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
	for(k in 1:nRepeat) {
		for(i in 1:nCptLen) {
			normalDataName = sapply(1:nPair, function(j) {return(paste(filename, "_Data_Normal_", j, "_", k, ".txt", sep=""))})
			tumorDataName = sapply(1:nPair, function(j) {return(paste(filename, "_Data_Tumor_Spike", cptLen[i], "_", j, "_", k, ".txt", sep=""))})
			dataFileName = c(normalDataName, tumorDataName)
			sampleName = c(rep(paste(filename, "_N","_", k, sep=""), nPair), rep(paste(filename, "_T", cptLen[i], "_", k, sep=""), nPair))
			sampleType = c(rep("Normal", nPair), rep("Tumor", nPair))
			metaMat = cbind(dataFileName, sampleName, sampleType)
			metaName = paste(filename, "_Meta_", cptLen[i], "_", k, ".txt", sep="")
			colnames(metaMat) = c("File", "Sample", "Type")
			write.table(metaMat, file=metaName, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
		}
		for(j in 1:nPair) {
			dataFileName = paste(filename, "_Data_Normal_", j, "_", k, ".txt", sep="")
			controlSim = nhppSimulate(controlRates)
			nLength = length(controlSim)
			strandSim = sample(0:1, nLength, replace=TRUE)
			chromosomeSim = rep(chromosomeN, nLength)
			dataMat = cbind(chromosomeSim, controlSim, strandSim)
			write.table(dataMat, file=dataFileName, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
		for(i in 1:nCptLen) {
			newSpikeRate = nhppSpikeConstWindow(controlRates, nSpike, cptLen[i], minGain, maxGain, minLoss, maxLoss, pGain)
			newSpikeMat = newSpikeRate$spikeMat
			spikeName = paste(filename, "_spikeMat_", cptLen[i], "_", k, ".txt", sep="")
			write.table(newSpikeMat, file=spikeName, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
			for(j in 1:nPair) {
				dataFileName = paste(filename, "_Data_Tumor_Spike", cptLen[i], "_", j, "_", k, ".txt", sep="")
				caseSim = nhppSimulate(newSpikeRate$caseRates)
				nLength = length(caseSim)
				strandSim = sample(0:1, nLength, replace=TRUE)
				chromosomeSim = rep(chromosomeN, nLength)
				dataMat = cbind(chromosomeSim, caseSim, strandSim)
				write.table(dataMat, file=dataFileName, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
				#newSpikeCBSN5 = ScanCBS(controlSim, newSpikeSim, statistic="normal", takeN=5, minStat=8, maxNCut=50, timing=TRUE)
				#newCBSRes = list(spikeCBSN5=newSpikeCBSN5, spikeMat=newSpikeRate$spikeMat)
				#simCBSRes[[counter]] = newCBSRes
			}
		}
	}
}

