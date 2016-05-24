# TODO: Add comment
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

testCreateExpressionDataset <- function(){
	
	cat("--- testCreateExpressionDataset: --- \n")
	
	stopifnot(all.equal(pData(eset),expDesign))
	stopifnot(all.equal(sampleNames(eset),colnames(m)))
	stopifnot(all.equal(nrow(exprs(eset)),nrow(m)))
	cat("--- testCreateExpressionDataset: PASS ALL TEST --- \n")
	
}


testGetAllEBayes <- function(){
	
	
	cat("--- testGetAllEBayes: --- \n")
	p <- getAllEBayes(eset)
	stopifnot( mean(p[,"A"]) > mean(p[,"B"]) )
	p <- getAllEBayes(eset)
	stopifnot( mean(p[,"A"]) > mean(p[,"B"]) )
	pAdj <- getAllEBayes(eset, adjust=T)
	
	stopifnot( mean(pAdj[,"A"]) > mean(p[,"A"]) )
	cat("--- testGetAllEBayes: PASS ALL TEST --- \n")
	
	print(head(p))
	
}

testGetRatios <- function(){
	
	cat("--- testGetRatios: --- \n")
	
	r <- getRatios(eset)
	
	exprs(eset)[1,c(5:6,1:2)]
	#  C_rep_1  C_rep_2  A_rep_1  A_rep_2 
	# 9.963852 9.966003 9.965568 9.967214 
	
	
	## NOTE: C is control
	stopifnot(all.equal(r[1,1],log2(median(exprs(eset)[1,1:2]) / median(exprs(eset)[1,5:6]))))
	
	stopifnot(mean(r[,"A"]) < mean(r[,"B"]))
	stopifnot(all.equal(r,getRatios(eset,method="mean")))
	r <- getRatios(eset)
	stopifnot(all.equal(mean(r[,"B"]),mean(r[,"B"])))
	
	cat("--- testGetRatios: PASS ALL TEST --- \n")
	
}

testGetAllCV <- function(){
	
	cat("--- testGetAllCV: --- \n")
	
	cv <- getAllCV(eset)
	
	stopifnot(all.equal(cv[1,"A"] , (sd(exprs(eset)[1,unlist(pData(eset)$condition == "A")]) / mean(exprs(eset)[1,unlist(pData(eset)$condition == "A")]))))
	stopifnot(all.equal(cv[200,"C"] , sd(exprs(eset)[200,unlist(pData(eset)$condition == "C")]) / mean(exprs(eset)[200,unlist(pData(eset)$condition == "C")])))
	
	cat("--- testGetAllCV: PASS ALL TEST --- \n")
}

testGlobalNormalize <- function(){
	
	cat("--- testGlobalNormalize: --- \n")
	
	globalNormFactors <- getGlobalNormFactors(eset)
	### add normalization factors to ExpressionSet
	pData(eset) <- cbind(pData(eset),globalNormFactors)
	esetNorm <- globalNormalize(eset,globalNormFactors)
	
	stopifnot(all.equal(as.vector(unlist(apply(exprs(esetNorm),2,sum))),as.vector(rev(unlist(apply(exprs(esetNorm),2,sum))))))
	
	stopifnot(pData(esetNorm)$normFactor[1] == 1)
	stopifnot(pData(esetNorm)$normFactor[2] != 1)
	
	cat("--- testGlobalNormalize: PASS ALL TEST --- \n")
	
	# Should generate error and stop	
#	fData(eset)$isNormAnchor <- rep(F,nrow(eset))
#	esetNorm <- normalizeIntensities(eset)
	
}

testGetSignalPerCondition <- function(){
	
	cat("--- testGetSignalPerCondition: --- \n")
	stopifnot(sum(getSignalPerCondition(eset,method="min")[,"A"] <= getSignalPerCondition(eset,method="max")[,"A"]) == nrow(eset))
	stopifnot(sum(getSignalPerCondition(eset,method="min")[,"C"] <= getSignalPerCondition(eset,method="median")[,"C"]) == nrow(eset))
	cat("--- testGetSignalPerCondition: PASS ALL TEST --- \n")
}

testBaselineIntensity <- function(){
	
	cat("--- testBaselineIntensity: --- \n")
	
	allInt <- as.vector(unlist(exprs(eset)))
	bl <- round(getBaselineIntensity(allInt,promille=5),2)
	#hist(allInt)
	#abline(v=bl,lwd=2)
	stopifnot(all.equal(bl[[1]] , 997.82) )
	
	
	set.seed(1234)
	suppressWarnings(allInt2 <- log10(rnorm(1000,1,1)))
	bl2 <- round(getBaselineIntensity(allInt2,promille=5),2)
	stopifnot(all.equal(-1.95,bl2[[1]]))
	#hist(allInt2)
	#abline(v=bl2,lwd=2)
	
	cat("--- testBaselineIntensity: PASS ALL TEST --- \n")
	
}

testRollUp <- function(){
	
	cat(" --- testRollUp --- \n")
	
	rollUpEset1 <- rollUp(eset,featureDataColumnName= c("proteinName"), method=c("sum"))
	stopifnot( length( unique( fData(eset)$proteinName ) ) == nrow(rollUpEset1)) 
	
	rollUpEset2 <- rollUp(eset ,featureDataColumnName= c("ptm"), method=c("sum"))
	stopifnot( length( unique( fData(eset)$ptm ) ) == nrow(rollUpEset2)) 
	
	rollUpEset3 <- rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("ptm"), method=c("mean"))
	stopifnot( length( unique( fData(eset)$ptm ) ) == nrow(rollUpEset3)) 
	
	print(exprs(rollUpEset2))
	
	stopifnot(all.equal(sum(exprs(rollUpEset2)),sum(exprs(rollUpEset1)))) ### test sum
	stopifnot(sum(exprs(rollUpEset1)) != sum(exprs(rollUpEset3))) ### test mean
	
	rollUpEset4 <- rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("proteinName"), method=c("top3"))
	stopifnot(sum(exprs(rollUpEset3)) != sum(exprs(rollUpEset4))  ) ### test top 3
	
	cat(" --- testRollUp: PASS ALL TEST  --- \n")
	
#	progenesisFeatureCsvFile3 <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/inst/testData/2014/peptides2.csv"
#	d <- parseProgenesisFeatureCsv(file=progenesisFeatureCsvFile3,expDesign=getExpDesignProgenesisCsv(progenesisFeatureCsvFile3))
#	system.time(e <- rollUp(d, method="sum", isProgressBar=T, featureDataColumnName= c("peptide")))
#	system.time(e <- rollUpDT(d, method="sum",  featureDataColumnName= c("peptide")))
#	
	esetTmp <- parseProgenesisPeptideMeasurementCsv(progenesisPeptideMeasurementCsvFile1,expDesign= getExpDesignProgenesisCsv(progenesisPeptideMeasurementCsvFile1))
	
	fData(esetTmp)
	
	rollUpEsetProteinAllAccessions <- rollUp(esetTmp,featureDataColumnName= c("proteinName"), method=c("sum"))
	stopifnot(fData(rollUpEsetProteinAllAccessions)$allAccessionsTMP[68] == "sp|E9PAV3|NACAM_HUMAN;sp|Q9BZK3|NACP1_HUMAN")

}

testTopX <- function(){
	
	cat(" --- testTopX --- \n")
	
	entryData1  <- data.frame(t(matrix(c(1,1,1,3,3,3,2,2,2,5,5,5),ncol=4)))
	rownames(entryData1) <- paste("peptide",1:nrow(entryData1),sep="_")
	#           X1 X2 X3
	# peptide_1  1  1  1
	# peptide_2  3  3  3
	# peptide_3  2  2  2
	# peptide_4  5  5  5
	stopifnot(all.equal(rep(10/3,3) , as.vector(unlist(getTopX(entryData1))) ))
	
	entryData2 <- data.frame(t(matrix(c(1,1,1,3,3,3,2,2,2,5,5,NA),ncol=4)))
	rownames(entryData2) <- paste("peptide",1:nrow(entryData2),sep="_")
	
	#           X1 X2 X3
	# peptide_1  1  1  1
	# peptide_2  3  3  3
	# peptide_3  2  2  2
	# peptide_4  5  5 NA
	stopifnot(all.equal(c(4,4,3) ,  as.vector(unlist(getTopX(entryData2,topX=2)))))
	
	# 1 row
	stopifnot(all.equal(rep(1,3),as.vector(unlist(getTopX(entryData1[1,])))))
	
	# 1 col
	stopifnot(all.equal(getTopX(entryData1)[[1]],getTopX(entryData1[,1])))
	
	top1 <- apply(exprs(rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("proteinName"), method=c("top1"))),1,sum)
	top3 <- apply(exprs(rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("proteinName"), method=c("top3"))),1,sum)
	meanInt <- apply(exprs(rollUp(eset[!fData(eset)$isFiltered,] ,featureDataColumnName= c("proteinName"), method=c("mean"))),1,sum)
	
	stopifnot(sum(top1 >=  top3 ) == length(top3))
	stopifnot(sum(top1) > sum(top3))
	stopifnot(sum(top1) > sum(meanInt))
	stopifnot(all.equal(sum(top3),sum(meanInt)))
	
	cat(" --- testTopX: PASS ALL TEST  --- \n")
	
}

testGetIBAQEset <- function(){
	
	cat(" --- testGetIBAQEset --- \n")
	
	# read protein fasta
	proteinDB <- read.fasta(fastaFile,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
	
	iBaqEset <- getIBAQEset(eset,proteinDB=proteinDB)
	stopifnot(round(exprs(iBaqEset))[1,1] == 125)
	stopifnot(round(exprs(iBaqEset))[2,2] == 38)
	
	cat(" --- testGetIBAQEset: PASS ALL TEST --- \n")
}


testGetLoocvFoldError <- function(){
	
	cat(" --- testGetLoocvFoldError --- \n")
	#plotCalibrationCurve(fit)
	stopifnot(  round(sum(getLoocvFoldError(absEstSimData))) == -8)
	cat(" --- testGetLoocvFoldError: PASS ALL TEST --- \n")
}

testSqNormalize <- function(){
	
	cat(" --- testSqNormalize --- \n")
	
	stopifnot(nrow(sqNormalize(eset, method = "global")) == 900)
	esetTmp <- parseProgenesisFeatureCsv(file=progenesisFeatureCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisFeatureCsvFile1))
	stopifnot(nrow(sqNormalize(esetTmp, method = "rt")) == 496)
	
	cat(" --- testSqNormalize: PASS ALL TEST --- \n")
}
	

testRtNormalize <- function(){
	

	cat(" --- testRTNormalize --- \n")
	esetTmp <- parseProgenesisFeatureCsv(file=progenesisFeatureCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisFeatureCsvFile1))
	rtNormFactors <- getRTNormFactors(esetTmp, minFeaturesPerBin=100)
	stopifnot(nrow(rtNormalize(esetTmp,rtNormFactors)) == 97)
	
	# Stop if rtNormFactors doesn't cover all retention times.
	#rtNormalize(esetTmp,rtNormFactors[1:2,])
	
	cat(" --- testRTNormalize: PASS ALL TEST --- \n")
}

testRemoveOutliers <- function(){
	
	cat(" --- testRemoveOutliers --- \n")
	set.seed(1234)
	stopifnot(sum(is.na(removeOutliers(c(-10,  rnorm(100), 10)))) == 2)
	cat(" --- testRemoveOutliers: PASS ALL TEST --- \n")
}

testPerFeatureNormalization <- function(){
	
	cat(" --- testPerFeatureNormalization: --- \n")
	
	normFactors <- exprs(eset)[1:10,1:3]
	colnames(normFactors) <- c("A","B","C")
	rownames(normFactors) <- fData(eset)[1:10,]$proteinName
	normFactors[is.finite(normFactors)] <- 1
	eNorm <-  perFeatureNormalization(eset,normFactors)
	
	stopifnot(all.equal(exprs(eNorm)[1,1] , (exprs(eset)[1,1]-1)))
	stopifnot(all.equal(exprs(eNorm)[1,2] , (exprs(eset)[1,2]-1)))
	stopifnot(all.equal(exprs(eNorm)[13,2] , (exprs(eset)[13,2])))
	
	normFactors[,2] <- 200
	eNorm <-  perFeatureNormalization(eset,normFactors)
	stopifnot(all.equal(exprs(eNorm)[1,1] , (exprs(eset)[1,1]-1)))
	stopifnot(all.equal(exprs(eNorm)[1,3] ,  (exprs(eset)[1,3]-200)))
	stopifnot(all.equal(exprs(eNorm)[1,4] ,  (exprs(eset)[1,4]-200)))
	stopifnot(all.equal(exprs(eNorm)[1,5] ,  (exprs(eset)[1,5]-1)))
	
	normFactors[,3] <- 0
	eNorm <-  perFeatureNormalization(eset,normFactors)
	stopifnot(all.equal(exprs(eNorm)[1,5] , exprs(eset)[1,5]))
	stopifnot(all.equal(exprs(eNorm)[2,6] , exprs(eset)[2,6]))
	
	normFactors[3,] <- 1000 
	eNorm <-  perFeatureNormalization(eset,normFactors)
	stopifnot(all.equal(exprs(eNorm)[3,] , (exprs(eset)[3,]-1000)))
	stopifnot(all.equal(exprs(eNorm)[20:50,] , exprs(eset)[20:50,]))
	
	# re-order normFactors columns
	eNorm <-  perFeatureNormalization(eset,normFactors[,rev(colnames(normFactors))])
	stopifnot(all.equal(exprs(eNorm)[3,] , (exprs(eset)[3,]-1000)))
	stopifnot(all.equal(exprs(eNorm)[20:50,] , exprs(eset)[20:50,]))
	
	cat(" --- testPerFeatureNormalization: PASS ALL TEST --- \n")
	
	#coveredPeptideSel <- fData(eset)$proteinName %in% rownames(normFactors)
	#exprs(eset)[coveredPeptideSel,]	<- exprs(eset)[coveredPeptideSel, ] - normFactors[as.character(fData(eset)[coveredPeptideSel,]$proteinName),pData(eset)$condition]
	
}


testStandardise <- function(){
	
	cat(" --- testStandardise: --- \n")
	
	d <- data.frame(rnorm(10000,5,20), rnorm(10000,10,10))
	dStd <- standardise(d[,1])
	
	stopifnot(round(mean(dStd)) == 0 )
	stopifnot(round(sd(dStd)) == 1 )
	
	dStd <- standardise(d)
	stopifnot(round(apply(dStd,2,mean)[1]) == 0)
	stopifnot(round(apply(dStd,2,sd)[1]) == 1)

	cat(" --- testStandardise: PASS ALL TEST --- \n")
	
}

testGetMaxIndex <-function(){
	
	cat(" --- testGetMaxIndex:  --- \n")
	d <- data.frame(s=c(NA,NA,NA,NA,1,1:4),lab=sort(rep(c("A","B","C"),3)))
	DT <- data.table(d)
	setkey(DT,lab)
	
	
	stopifnot(all.equal(c(1,5,9) , DT[, .I[getMaxIndex(s)], by=lab ]$V1  ))
	cat(" --- testGetMaxIndex: PASS ALL TEST  --- \n")
	
}

### TEST FUNCTIONS END

### TESTS
testCreateExpressionDataset()
testGetAllEBayes()
testGetRatios()
testGetAllCV()
testGlobalNormalize()
#testNormalise()
testRtNormalize()
testGetSignalPerCondition()
testBaselineIntensity()
testRollUp()
testTopX()

testGetLoocvFoldError()

testRemoveOutliers()
testPerFeatureNormalization()
testStandardise()
testGetMaxIndex()

testGetIBAQEset()

