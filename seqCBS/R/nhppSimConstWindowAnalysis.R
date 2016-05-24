nhppSimConstWindowAnalysis <-
function(filePrefix, chromosomeN, distMetric=c(20,50,100,150,200,300,500,1000), cptLen=c(3,5,8,12,15,20,30,50,100), nPair=2, nRepeat=10, statistic="normal", grid.size="auto", takeN=5, maxNCut=60, minStat=5, verbose=FALSE, timing=TRUE, hasRun=FALSE, width=12, height=6) {
	nCptLen = length(cptLen)
	tauHatCBS = vector("list", nCptLen)
	tauHatSeg = vector("list", nCptLen)
	simCBS = vector("list", nCptLen)
	trueTau = vector("list", nCptLen)
	tauCBSDist = vector("list", nCptLen)
	tauSegDist = vector("list", nCptLen)
	CBSMatchDist = vector("list", nCptLen)
	SegMatchDist = vector("list", nCptLen)
	CBSRecall = matrix(0, nrow=length(cptLen), ncol=length(distMetric))
	SegRecall = matrix(0, nrow=length(cptLen), ncol=length(distMetric))
	CBSPrecision = matrix(0, nrow=length(cptLen), ncol=length(distMetric))
	SegPrecision = matrix(0, nrow=length(cptLen), ncol=length(distMetric))
	CBSFMeasure = matrix(0, nrow=length(cptLen), ncol=length(distMetric))
	SegFMeasure = matrix(0, nrow=length(cptLen), ncol=length(distMetric))
	trueTauMeanSigLen = vector("list", nCptLen)
	nTrueTau = rep(0, length(cptLen))
	nCBSCall = rep(0, length(cptLen))
	nSegCall = rep(0, length(cptLen))
	CBSTime = rep(0, length(cptLen))
	for(i in 1:nCptLen) {
		tauHatCBS[[i]] = vector("list", nRepeat)
		tauHatSeg[[i]] = vector("list", nRepeat)
		simCBS[[i]] = vector("list", nRepeat)
		trueTau[[i]] = vector("list", nRepeat)
		tauCBSDist[[i]] = vector("list", nRepeat)
		tauSegDist[[i]] = vector("list", nRepeat)
		CBSMatchDist[[i]] = vector("list", nRepeat)
		SegMatchDist[[i]] = vector("list", nRepeat)
	}
	for(i in 1:nCptLen) {
		for(j in 1:nRepeat) {
			metaName = paste(filePrefix, "_Meta_", cptLen[i], "_", j, ".txt", sep="")
			spikeName = paste(filePrefix, "_spikeMat_", cptLen[i], "_", j, ".txt", sep="")
			CBSResName = paste(filePrefix, "_CBS_", cptLen[i], "_", j, "_", statistic, ".Rdata", sep="")
			trueSpike = as.matrix(read.delim(spikeName, header=FALSE))
			seqs = readInput(metaName, formatName="Chiang")
			normalSeq = seqs$normalSeq$seqF
			tumorSeq = seqs$tumorSeq$seqF
			rm(seqs)
			combL = CombineCaseControlC(tumorSeq, normalSeq)$combL
			
			trueTau2 = floor(sort(trueSpike[,c(3,4)]))
			trueTau[[i]][[j]] = findInterval(trueTau2, combL)
			trueSpikeL = t(apply(trueSpike[,c(3,4)], 1, findInterval, vec=combL))
			trueTauMaxSigLen = max(trueSpikeL[,2]-trueSpikeL[,1])
			trueTauSigLen2 = diff(trueTau[[i]][[j]])
			trueTauSigLen = trueTauSigLen2[trueTauSigLen2<=trueTauMaxSigLen]
			trueTauMeanSigLen[[i]] = c(trueTauMeanSigLen[[i]], trueTauSigLen)

			if(!hasRun) {
				simCBS[[i]][[j]] = ScanCBS(tumorSeq, normalSeq, statistic, grid.size, takeN, maxNCut, minStat, verbose, timing)
				simCBSRes=simCBS[[i]][[j]]
				save(simCBSRes, file=CBSResName)
			}
			else {
				load(CBSResName)
				simCBS[[i]][[j]] = simCBSRes
			}
			ScanCBSSimPlot(tumorSeq, normalSeq, simCBS[[i]][[j]], trueTau2, trueSpike, filename=paste(filePrefix, "_CBSPlot_",  cptLen[i], "_", j, "_", statistic, sep=""), mainTitle=paste(filePrefix, "_CBSPlot_",  cptLen[i], "_", j, "_", statistic, sep=""))
			tauHatCBS[[i]][[j]] = sort(simCBS[[i]][[j]]$statHat[,3:4])
			tauCBSDist[[i]][[j]] = t(sapply(trueTau[[i]][[j]], function(x) {sapply(tauHatCBS[[i]][[j]], function(y){return(as.numeric(abs(x-y)))})}))
			segName = paste("JSSim", cptLen[i], "_", j, "_seg_W_400_initFP_1000_finalFP_10.txt", sep="")
			tauHatSeg[[i]][[j]] = findInterval(SegSeqResProcess(segName)[[1]], combL)
			tauSegDist[[i]][[j]] = t(sapply(trueTau[[i]][[j]], function(x) {sapply(tauHatSeg[[i]][[j]], function(y){return(as.numeric(abs(x-y)))})}))
			nTrueTau3 = nrow(tauSegDist[[i]][[j]])
			nCBSCall3 = ncol(tauCBSDist[[i]][[j]])
			nSegCall3 = ncol(tauSegDist[[i]][[j]])
			if(nCBSCall3 < nTrueTau3)	tauCBSDist[[i]][[j]] = t(tauCBSDist[[i]][[j]])
			if(nSegCall3 < nTrueTau3)	tauSegDist[[i]][[j]] = t(tauSegDist[[i]][[j]])
			CBSMatch = solve_LSAP(tauCBSDist[[i]][[j]])
			CBSMatch = cbind(seq_along(CBSMatch), CBSMatch)
			CBSMatchDist[[i]][[j]] = apply(CBSMatch, 1, function(x) {return(tauCBSDist[[i]][[j]][x[1], x[2]])})
			SegMatch = solve_LSAP(tauSegDist[[i]][[j]])
			SegMatch = cbind(seq_along(SegMatch), SegMatch)
			SegMatchDist[[i]][[j]] = apply(SegMatch, 1, function(x) {return(tauSegDist[[i]][[j]][x[1], x[2]])})
			if(verbose) {
				print(c(nTrueTau2, nCBSCall2, nSegCall2))
				print(CBSMatchDist[[i]][[j]])
				print(SegMatchDist[[i]][[j]])
			}
		}
		for(j in 1:nRepeat) {
			nTrueTau2 = length(trueTau[[i]][[j]])
			nCBSCall2 = length(tauHatCBS[[i]][[j]])
			nSegCall2 = length(tauHatSeg[[i]][[j]])
			nTrueTau[i] = nTrueTau[i]+nTrueTau2
			nCBSCall[i] = nCBSCall[i]+nCBSCall2
			nSegCall[i] = nSegCall[i]+nSegCall2
			nCBSTruePos = sapply(distMetric, function(x) {sum(CBSMatchDist[[i]][[j]] <= x)})
			nSegTruePos = sapply(distMetric, function(x) {sum(SegMatchDist[[i]][[j]] <= x)})
			CBSRecall[i,] = CBSRecall[i,]+nCBSTruePos
			CBSPrecision[i,] = CBSPrecision[i,]+nCBSTruePos
			SegRecall[i,] = SegRecall[i,]+nSegTruePos
			SegPrecision[i,] = SegPrecision[i,]+nSegTruePos
			CBSTime[i] = CBSTime[i] + simCBS[[i]][[j]]$timingRes$timeCBSTotal/nRepeat
		}
		CBSRecall[i,] = CBSRecall[i,]/nTrueTau[i]
		CBSPrecision[i,] = CBSPrecision[i,]/nCBSCall[i]
		SegRecall[i,] = SegRecall[i,]/nTrueTau[i]
		SegPrecision[i,] = SegPrecision[i,]/nSegCall[i]
		nTrueTau[i] = nTrueTau[i]/nRepeat
		nCBSCall[i] = nCBSCall[i]/nRepeat
		nSegCall[i] = nSegCall[i]/nRepeat
	}
	CBSFMeasure = 2*CBSPrecision*CBSRecall/(CBSPrecision+CBSRecall)
	SegFMeasure = 2*SegPrecision*SegRecall/(SegPrecision+SegRecall)
	trueTauMeanSigLen = sapply(trueTauMeanSigLen, mean, trim=0.05)
	
	nDistMetric = length(distMetric)
	plotFileName = paste(filePrefix, "_SampleCalling_", statistic, ".pdf", sep="")
	pdf(plotFileName, width=width, height=height)
	par(mfrow=c(1,2))
	plot(x=cptLen*20, y=trueTauMeanSigLen, type="p", pch=20, cex=2, xlab="Simulation Signal Length", ylab="Mean Distance Between True Signal Change Points", main="Distance between True Signal Change Points")
	matplot(x=log(cbind(trueTauMeanSigLen, trueTauMeanSigLen), base=10), y=cbind(nCBSCall, nSegCall), type="p", pch=c(20,20), cex=2, xlab="Log Mean Distance Between True Signal Change Points", ylab="Number of Cpts Called", main="Number of Change Points Called")
	legend("bottomright", legend=c("SeqCBS", "SegSeq"), pch=c(20,20), col=1:2)
	dev.off()

	for(k in 1:nDistMetric) {
		plotFileName1 = paste(filePrefix, "_Performance1_", statistic, "_Metric", distMetric[k], ".pdf", sep="")
		pdf(plotFileName1, width=width, height=height)
		par(mfcol=c(1,2))
		matplot(x=log(cbind(trueTauMeanSigLen,trueTauMeanSigLen),base=10), y=cbind(CBSRecall[,k], SegRecall[,k]), type="p", pch=c(20,20), cex=2, xlab="Log Mean Signal Length", ylab="Recall Rate", ylim=0:1, main=paste("Recall Rate, Metric=",distMetric[k],sep=""))
		abline(h=c(mean(CBSRecall[,k]), mean(SegRecall[,k])))
		legend("bottomright", legend=c("SeqCBS", "SegSeq"), pch=c(20,20), col=1:2)
		matplot(x=log(cbind(trueTauMeanSigLen,trueTauMeanSigLen),base=10), y=cbind(CBSPrecision[,k], SegPrecision[,k]), type="p", pch=c(20,20), cex=2, xlab="Log Mean Signal Length", ylab="Precision Rate", ylim=0:1, main=paste("Precision Rate, Metric=",distMetric[k],sep=""))
		abline(h=c(mean(CBSPrecision[,k]), mean(SegPrecision[,k])))
		legend("bottomright", legend=c("SeqCBS", "SegSeq"), pch=c(20,20), col=1:2)
		dev.off()
		
		plotFileName2 = paste(filePrefix, "_Performance2_", statistic, "_Metric", distMetric[k], ".pdf", sep="")
		pdf(plotFileName2, width=width, height=height)
		par(mfcol=c(1,2))		
		matplot(x=cbind(CBSPrecision[,k], SegPrecision[,k]), y=cbind(CBSRecall[,k], SegRecall[,k]), type="p", pch=c(20,20), cex=2, xlab="Precision", ylab="Recall Rate", xlim=0:1, ylim=0:1, main=paste("Precision vs Recall, Metric=",distMetric[k],sep=""))
		legend("bottomright", legend=c("SeqCBS", "SegSeq"), pch=c(20,20), col=1:2)
		matplot(x=log(cbind(trueTauMeanSigLen,trueTauMeanSigLen),base=10), y=cbind(CBSFMeasure[,k], SegFMeasure[,k]), type="p", pch=c(20,20), cex=2, xlab="Log Mean Signal Length", ylab="F Measure", ylim=0:1, main=paste("F Measure, Metric=",distMetric[k],sep=""))
		abline(h=c(mean(CBSFMeasure[,k]), mean(SegFMeasure[,k])))
		legend("bottomright", legend=c("SeqCBS", "SegSeq"), pch=c(20,20), col=1:2)
		dev.off()		
	}
	return(list(simCBS=simCBS, CBSMatchDist=CBSMatchDist, SegMatchDist=SegMatchDist, CBSRecall=CBSRecall, SegRecall=SegRecall, CBSPrecision=CBSPrecision, SegPrecision=SegPrecision, CBSFMeasure=CBSFMeasure, SegFMeasure=SegFMeasure, trueTauMeanSigLen=trueTauMeanSigLen, nTrueTau=nTrueTau, nCBSCall=nCBSCall, nSegCall=nSegCall, CBSTime=CBSTime))
}
