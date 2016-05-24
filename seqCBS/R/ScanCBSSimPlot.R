ScanCBSSimPlot <-
function(cases, controls, CBSObj, trueTau, SpikeMat, filename, mainTitle, CIObj=NULL, length.out=10000, localWindow=0.5*10^5, localSeparatePlot=TRUE, smoothF=0.025, xlabScale=10^6, width=12, height=18) {
	p = length(cases)/(length(cases)+length(controls))
	maxCase = max(cases)
	maxControl = max(controls)
	maxVal = max(c(maxCase, maxControl))
	cpts = matrix(CBSObj$statHat[,c(1,2,5)], nrow=nrow(CBSObj$statHat))
	tauHatFull = CBSObj$tauHat
	tauHat = tauHatFull[c(-1, -length(tauHatFull))]
	relCN = CBSObj$relCN
	relCN[relCN <= 0] = 1
	ylims = c(min(c(0, cpts[,3])), max(c(0, cpts[,3])))
	if(!is.null(CIObj)) {
		CIBounds = CIObj$CIRes[3:4,]
		CIL = as.numeric(CIObj$CIRes[5,])
		CIU = as.numeric(CIObj$CIRes[6,])
		relCNCIL = CIL/(1-CIL+10^(-5))/(p/(1-p))
		relCNCIU = CIU/(1-CIU+10^(-5))/(p/(1-p))
	}
	
	grid.fix = seq(1, maxVal, length.out=length.out)
	gridSize = grid.fix[2]-grid.fix[1]
	casesCountInGrid = getCountsInWindow(cases, 0, maxVal, gridSize, sorted=FALSE)
	casesCountInGridSmooth = lowess(x=grid.fix, y=casesCountInGrid, smoothF)
	casesCountInGridSmooth$y[casesCountInGridSmooth$y<0] = 0
	controlCountInGrid = getCountsInWindow(controls, 0, maxVal, gridSize, sorted=FALSE)
	controlCountInGridSmooth = lowess(x=grid.fix, y=controlCountInGrid, smoothF)
	controlCountInGridSmooth$y[controlCountInGridSmooth$y<0] = 0
	PInGrid = casesCountInGrid/(casesCountInGrid+controlCountInGrid)
	PInGrid[is.nan(PInGrid)]=0
	PInGridSmooth = lowess(x=grid.fix, y=PInGrid, smoothF)
	relCNInGrid = PInGrid/(1-PInGrid)/(p/(1-p))
	relCNInGrid[is.nan(relCNInGrid) | !is.finite(relCNInGrid) | relCNInGrid <= 0]=1
	relCNInGrid = log(relCNInGrid, base=2)
	relCNlims = c(min(min(log(relCN, base=2)), min(relCNInGrid)), max(max(log(relCN, base=2)), max(relCNInGrid)))
	tauHatInGrid = grid.fix[tauHat %/% gridSize]/xlabScale
	trueTauInGrid = grid.fix[trueTau %/% gridSize]/xlabScale
	gridYLims = c(min(log(casesCountInGrid+1) - log(controlCountInGrid+1)), log(max(casesCountInGrid, controlCountInGrid)))
	
	plotTauHatInd = c(min(min(cases),min(controls)), tauHat, maxVal) %/% gridSize
	plotTauHatInd = sapply(plotTauHatInd, function(x) {max(x,1)})
	plotTauHatInd = sapply(plotTauHatInd, function(x) {min(x,max(grid.fix))})
	plotTauHat = grid.fix[plotTauHatInd]/xlabScale
	
	## 1. Plot the chromosome global view
	pdf(paste(filename, ".pdf", sep=""), width=width, height=height)
	par(mfrow=c(3,1))
	
	plot(x=grid.fix/xlabScale, y=rep(0, length(grid.fix)), type="n", ylim=ylims, main=mainTitle, ylab="Statistic", xlab=paste("Base Pairs", xlabScale))
	for(i in 1:nrow(cpts)) {
		plotX = c(grid.fix[max(floor(cpts[i,1]/gridSize), 1)]/xlabScale, grid.fix[ceiling(cpts[i,2]/gridSize)]/xlabScale)
		lines(x=plotX, y=rep(cpts[i,3],2), lwd=3)
	}
	for(i in 1:nrow(SpikeMat)) {
		plotX = c(grid.fix[max(floor(SpikeMat[i,3]/gridSize), 1)]/xlabScale, grid.fix[ceiling(SpikeMat[i,4]/gridSize)]/xlabScale)
		lines(x=plotX, y=rep(SpikeMat[i,5],2), lwd=2, col=2)
	}
	abline(v=tauHatInGrid, lty=3, col=4)
	abline(v=trueTauInGrid, lty=3, col=2)
	
	matplot(x=grid.fix/xlabScale, y=log(cbind(casesCountInGridSmooth$y, controlCountInGridSmooth$y)+1), type="l", lty=c(1,1), col=c(2,1), main="Log Read Intensity", ylab="Read Intensity", xlab=paste("Base Pairs", xlabScale), ylim=gridYLims)
	points(x=grid.fix/xlabScale, y=log(casesCountInGrid+1) - log(controlCountInGrid+1), pch=".", col=1)
	abline(v=tauHatInGrid, lty=3, col=4)
	abline(v=trueTauInGrid, lty=3, col=2)
	legend("topright", c("case","control", "case-control"), pch=".", lty=c(1,1,0), col=c(2,1,1))
	
	plot(x=grid.fix/xlabScale, y=relCNInGrid, type="p", pch=20, ylim=relCNlims, main="Log Relative Copy Number", ylab="Log2 Relative CN", xlab=paste("Base Pairs", xlabScale))
	lines(x=plotTauHat, y=log(c(relCN, relCN[length(relCN)]), base=2), type="s", col="red")
	abline(v=tauHatInGrid, lty=3, col=4)
	abline(v=trueTauInGrid, lty=3, col=2)
	dev.off()
	
	## 2. Plot Local View of each Change Point
	nTauHat = length(tauHat)
	if(localSeparatePlot == FALSE) {
		nPlotCol = as.integer(sqrt(nTauHat/(height/width)))
		nPlotRow = ceiling(nTauHat/nPlotCol)
		pdf(paste(filename, "_localDetails.pdf", sep=""), width=width*2, height=height*2)
		par(mfrow=c(nPlotRow, nPlotCol))
	}
	for(i in 1:nTauHat) {
		if(localSeparatePlot) {
			pdf(paste(filename, "_local_", i, "_", tauHat[i], ".pdf", sep=""), width=width, height=height/2)
		}
		lBound = max(0, tauHat[i]-localWindow)
		rBound = min(maxVal, tauHat[i]+localWindow)
		localCas = cases[cases >= lBound & cases < rBound]
		localCon = controls[controls >= lBound & controls < rBound]
		grid.fix = seq(lBound, rBound, length.out=length.out/100)
		gridSize = grid.fix[2]-grid.fix[1]
		grid.mpt = grid.fix + gridSize/2
		CasCountInGrid = getCountsInWindow(localCas, lBound, rBound, gridSize, sorted=FALSE)
		ConCountInGrid = getCountsInWindow(localCon, lBound, rBound, gridSize, sorted=FALSE)
		pInGrid = CasCountInGrid/(CasCountInGrid+ConCountInGrid)
		pInGrid[is.nan(pInGrid)] = 0.0
		
		combLocalCasCon = CombineCaseControlC(localCas, localCon)
		plotReadRangeInd = combLocalCasCon$combL >= lBound & combLocalCasCon$combL <= rBound
		plotReadX = combLocalCasCon$combL[plotReadRangeInd]
		plotReadY = combLocalCasCon$combZ[plotReadRangeInd] > 0
		
		plotPX = cbind(tauHatFull[-length(tauHatFull)], tauHatFull[-1])
		pSegment = relCN*p/(1-p)/(1+relCN*p/(1-p))
		plotPY = cbind(pSegment, pSegment)
		
		if(!is.null(CIObj)) {
			localCIBounds = (CIBounds[1,] <= rBound) & (CIBounds[2,] >= lBound)
			localYLims = c(min(CIL[localCIBounds]), max(CIU[localCIBounds])) * c(0.8, 1.2)
			if(is.nan(localYLims[1]) || !is.finite(localYLims[1]))	localYLims[1] = 0
			if(is.nan(localYLims[2]) || !is.finite(localYLims[2]))	localYLims[2] = 1
		}
		else {
			localYLims = c(0,1)
		}
		
		plot(x=1, y=1, type="n", xlim=c(lBound, rBound), xaxt="n", ylim=localYLims, main=paste("Reads and Inference around", tauHat[i]), xlab="Base Pair Locations", ylab="p(case read)", cex.main=0.75, cex.lab=0.75, cex.axis=0.75)
		axis(side=1, at=plotReadX[!plotReadY], labels=FALSE, tcl=0.3)
		axis(side=3, at=plotReadX[plotReadY], labels=FALSE, tcl=0.3)
		axis(side=1, xaxp=c(lBound, rBound, 10), tcl=-0.5, cex.axis=0.75)
		if(is.null(CIObj)) {
			if(length(grid.mpt) != length(pInGrid)) {
				grid.mpt = grid.mpt[1:max(length(grid.mpt), length(pInGrid))]
				pInGrid = pInGrid[1:max(length(grid.mpt), length(pInGrid))]
			}
			points(x=grid.mpt, y=pInGrid, pch=20, col=3)
		}
		else {
			for(j in 1:ncol(CIBounds)) {
				lines(x=CIBounds[,j], y=rep(CIL[j],2), col="#AAAAAA", lwd=2)
				lines(x=CIBounds[,j], y=rep(CIU[j],2), col="#AAAAAA", lwd=2)
			}
		}
		for(j in 1:nrow(plotPY)) {
			lines(x=plotPX[j,], y=plotPY[j,], lwd=3)
		}
		abline(v=tauHat, lty=3, lwd=2, col="#AAAAAA")
		abline(v=trueTau, lty=2, lwd=2, col="#777777")
		if(localSeparatePlot) {
			dev.off()
		}
	}
	if(localSeparatePlot == FALSE)	dev.off()
}
