# Grapihcs.R unit tests
# 
# Author: ahrnee-adm
###############################################################################

### INIT
if(!grepl("SafeQuant\\.Rcheck",getwd())){
	setwd(dirname(sys.frame(1)$ofile))
}
source("initTestSession.R")
### INIT END


### TEST FUNCTIONS

test.getConditionColors <- function(){
	
	cat("--- test.getConditionColors: --- \n")
	stopifnot(nrow(.getConditionColors(eset)) == 3)
	stopifnot(length(.getConditionColors(eset)[pData(eset)$condition,]) == 6)
	cat("--- test.getConditionColors: PASS ALL TEST --- \n")
}

### TEST FUNCTIONS END

### TESTS

test.getConditionColors()

testPlotNbIdentificationsVsRT <- function(){
	
	file <- progenesisPeptideMeasurementFile1
	#file <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/inst/testData/QI_2.0/peptide_measurements5_MRuegg.csv"
	
	cat("--- testPlotNbIdentificationsVsRT: --- \n")
	esetTmp <- parseProgenesisPeptideMeasurementCsv(file,expDesign= getExpDesignProgenesisCsv(file))
	plotNbIdentificationsVsRT(esetTmp)
	cat("--- testPlotNbIdentificationsVsRT: PASS ALL TEST --- \n")
	
}

if(T){
	
	tmpFile <- paste(tempdir(),"/tmp.pdf",collapse="",sep="")
	
	pdf(tmpFile)
	# plots
	plotExpDesign(eset)
	
	# id plots 
	isDec <- isDecoy(fData(eset)$proteinName)
	qvals <- getIdLevelQvals(fData(eset)$idScore, isDec)
	idScore <- fData(eset)$idScore
	pMassError <- fData(eset)$pMassError
	
	plotScoreDistrib(idScore[!isDec]
			,idScore[isDec], main="plotScoreDistrib")
	plotIdScoreVsFDR(idScore,qvals,lwd=2, main="plotIdScoreVsFDR")
	plotROC(qvals,breaks=seq(0,0.2,length=50)
			,main=paste("plotROC Nb. Valid identificaitons: ",sum(qvals < 0.01),"\n( FDR ",0.01,")"))
	
	# precursor mass error
	plotPrecMassErrorDistrib(eset,pMassTolWindow=c(-10,10),main="plotPrecMassErrorDistrib")
	plotPrecMassErrorVsScore(eset, pMassTolWindow=c(-0.5,0.5),main="plotPrecMassErrorVsScore")
	
	# quant QC plots
	pairsAnnot(exprs(eset), main="pairsAnnot")
	pairsAnnot(getSignalPerCondition(eset), main="pairsAnnot")
	
	
	
	plotMSSignalDistributions(exprs(eset),col=COLORS, cex=1, cex.axis=1.5, cex.lab=1.5, ylab="binned count", xlab="AUC", main="plotMSSignalDistributions")
	barplotMSSignal(eset,cex.lab=1.5,main="barplotMSSignal")
	
	##quant differential abundance related plots
	
	### plot volcanos for all case control comparisons
	plotVolcano(sqa
			,main= "plotVolcano created from safeQuantAnalysis object"
			,cex.axis=1.2
			,cex.lab=1.2
			,adjust=F
	)
	
	plotVolcano(sqa
			,main= "plotVolcano created from safeQuantAnalysis object"
			,cex.axis=1.2
			,cex.lab=1.2
			,adjust=T
	)
	
	plotVolcano(data.frame( ratio=as.vector(unlist(sqa$ratio["B"]))
					,qValue=as.vector(unlist(sqa$qValue["B"]))
					,cv=apply(sqa$cv[c("B","C")],1,max,na.rm=T))
			,caseCondition="cond B"
			,controlCondition="cond C"
			,main= "plotVolcano created from data.frame"
			,cex.axis=1.2
			,cex.lab=1.2
	)
	
	hClustHeatMap(eset)
	hClustHeatMap(eset, dendogram="both")
	
	par(mfrow=c(2,2))
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=T,isAdjusted=T,main="plotNbValidDeFeaturesPerFDR UP")
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=F,isAdjusted=T,main="plotNbValidDeFeaturesPerFDR DOWN")
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=T,isAdjusted=F,main="plotNbValidDeFeaturesPerFDR UP")
	plotNbValidDeFeaturesPerFDR(sqa,upRegulated=F,log2RatioCufOff=log2(1),pvalRange=c(0,0.3),pvalCutOff=1, isLegend=F,isAdjusted=F,main="plotNbValidDeFeaturesPerFDR DOWN")
	par(mfrow=c(1,1))
	
	plotXYDensity(exprs(eset)[,1],exprs(eset)[,2], main="plotXYDensity")
	plotAbsEstCalibrationCurve(absEstSimDataFit, predictorName="log10(iBAQ)", main="plotCalibrationCurve")

	esetTmp <- parseProgenesisFeatureCsv(file=progenesisFeatureCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisFeatureCsvFile1))
	rtNormFactors <- getRTNormFactors(esetTmp, minFeaturesPerBin=100)
	plotRTNormSummary(esetTmp, main="plotRTNormSummary")
	plotRTNorm(rtNormFactors,esetTmp, samples=2, main="plotRTNorm")
	
	missinValueBarplot(eset)
	
	plotQValueVsPValue(sqa, lim=c(0,0.5))
		
	.correlationPlot(exprs(eset))
	
	graphics.off()
	
	

	cat("CREATED FILE", tmpFile, "\n")
	
}
### TESTS END





