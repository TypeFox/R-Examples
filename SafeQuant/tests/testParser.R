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

################## TMT ##################

testGetSkipLineNb <- function(){
	cat(" --- testGetSkipLineNb --- \n")
	stopifnot(40 ==  .getSkipLineNb(scaffoldTmt6PlexRawTestFile))
	stopifnot(40 ==  .getSkipLineNb(scaffoldTmt10PlexRawTestFile))
	cat(" --- testGetSkipLineNb: PASS ALL TEST --- \n")	
}

testGetNbPlex <- function(){
	cat(" --- testGetNbPlex --- \n")
	stopifnot(10 ==.getNbPlex(scaffoldTmt10PlexRawTestFile))
	stopifnot(6 ==.getNbPlex(scaffoldTmt6PlexRawTestFile))
	cat(" --- testGetNbPlex: PASS ALL TEST --- \n")	
}



################## LFQ ##################

testGetProgenesisCsvExpressionColIndices <- function(){
	
	cat(" --- testGetProgenesisCsvProteinIntColIndices --- \n")
	
	stopifnot(sum(14:17 %in% .getProgenesisCsvExpressionColIndices(progenesisProteinCsvFile1) ) == 4 )
	
	stopifnot(sum(31:48 %in% .getProgenesisCsvExpressionColIndices(progenesisFeatureCsvFile1) ) == 18 )
	
	stopifnot(.getProgenesisCsvExpressionColIndices(progenesisPeptideMeasurementCsvFile1)  == 18 )
	
	stopifnot(.getProgenesisCsvExpressionColIndices(progenesisPeptideMeasurementFractionatedCsvFile1)[1] == 21)
	
	cat(" --- testGetProgenesisCsvProteinIntColIndices: PASS ALL TEST --- \n")	
	
}

testGetExpDesignProgenesisCsv <- function(){
	
	cat(" --- testGetExpDesignProgenesisCsv --- \n")
	
	stopifnot(4 == nrow( getExpDesignProgenesisCsv(progenesisProteinCsvFile1) )) 
	
	stopifnot(18 == nrow(getExpDesignProgenesisCsv(progenesisFeatureCsvFile1) ))
	
	cat(" --- testGetExpDesignProgenesisCsv: PASS ALL TEST  --- \n")
	
}


testParseProgenesisProteinCsv <- function(){
	
	cat(" --- testParseProgenesisProteinCsv --- \n")
	
	eset <- parseProgenesisProteinCsv(file=progenesisProteinCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisProteinCsvFile1))
	stopifnot(983 == nrow(eset))
	stopifnot(6 == ncol(fData(eset)))
	stopifnot(4 == nrow(pData(eset)))
	
	expDesign <- getExpDesignProgenesisCsv(progenesisProteinCsvFile1)[c(3,2),]
	eset2 <- parseProgenesisProteinCsv(file=progenesisProteinCsvFile1,expDesign=expDesign)
	
	stopifnot(colnames(eset2) == row.names(pData(eset))[c(3,2)]   )
	### change expDesign
	
	cat(" --- testParseProgenesisProteinCsv: PASS ALL TEST  --- \n")
}

testParseProgenesisFeatureCsv <- function(){
	
	cat(" --- testParseProgenesisFeatureCsv --- \n")
	
	eset <- parseProgenesisFeatureCsv(file=progenesisFeatureCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisFeatureCsvFile1))
	stopifnot(97 == nrow(eset))
	#names(fData(eset))
	stopifnot(13 == ncol(fData(eset)))
	stopifnot(18 == nrow(pData(eset)))
	
#	expDesign <- getExpDesignProgenesisCsv(progenesisFeatureCsvFile1)[c(3,4:8),]
#	eset2 <- parseProgenesisPeptideCsv(file=progenesisFeatureCsvFile1,expDesign=expDesign)
#	
#	pData(eset2)
#	pData(eset)
	
	cat(" --- testParseProgenesisFeatureCsv: PASS ALL TEST  --- \n")
}

testParseScaffoldRawFile <- function(){
	
	
	cat(" --- testParseScaffoldRawFile --- \n")
	
	# 6 plex	
	
	expDesignTMTSixPlex <- data.frame(condition=sort(rep(c(1,2),3)),isControl=sort(rep(c(T,F),3),decreasing=T) )
	eset <- parseScaffoldRawFile(scaffoldTmt6PlexRawTestFile,expDesign=expDesignTMTSixPlex)
	stopifnot(329 ==  nrow(eset))
	stopifnot("sp|P38711|RS27B_YEAST,sp|P35997|RS27_YEAST" ==  fData(eset)$proteinName[1])
	stopifnot("sp|P38711|RS27B_YEAST" == fData(parseScaffoldRawFile(scaffoldTmt6PlexRawTestFile,expDesign=expDesignTMTSixPlex,keepFirstAcOnly=T))$proteinName[1])
	
	expDesignTMTSixPlex <- data.frame(condition=sort(rep(c(1,2),3)),isControl=sort(rep(c(T,F),3),decreasing=T) )
	rownames(expDesignTMTSixPlex) <- rev(1:6)
	eset2 <- parseScaffoldRawFile(scaffoldTmt6PlexRawTestFile,expDesign=expDesignTMTSixPlex[1:5,])
	stopifnot(all.equal(exprs(eset)[1,6], exprs(eset2)[1,1]))
	
	# 10 plex	
	expDesignTMTTenPlex <- data.frame(condition=sort(rep(c(1:5),2)),isControl=c(T,T,rep(F,8)) )
	eset3 <- parseScaffoldRawFile(scaffoldTmt10PlexRawTestFile,expDesign=expDesignTMTTenPlex, keepFirstAcOnly=T)
	stopifnot(958 ==  nrow(eset3))
	cat(" --- testParseScaffoldRawFile: PASS ALL TEST --- \n")
	
}

testGetFileType <- function(){
	
	cat(" --- testGetFileType --- \n")
	
	stopifnot(.getFileType(scaffoldTmt6PlexRawTestFile) == "ScaffoldTMT")
	stopifnot(.getFileType(scaffoldTmt10PlexRawTestFile) == "ScaffoldTMT")
	stopifnot(.getFileType(progenesisFeatureCsvFile1) == "ProgenesisFeature")
	stopifnot(.getFileType(progenesisProteinCsvFile1) == "ProgenesisProtein")
	stopifnot(.getFileType(progenesisPeptideMeasurementCsvFile1) == "ProgenesisPeptide")
	stopifnot(.getFileType(maxQuantProteinFileTxt) == "MaxQuantProteinGroup")

	cat(" --- testGetFileType: PASS ALL TEST --- \n")
	
}

testParseProgenesisPeptideMeasurementCsv <- function(){
	
	cat(" --- testParseProgenesisPeptideMeasurementCsv --- \n")
	eset <- parseProgenesisPeptideMeasurementCsv(progenesisPeptideMeasurementCsvFile1,expDesign= getExpDesignProgenesisCsv(progenesisPeptideMeasurementCsvFile1))
	stopifnot(ncol(exprs(eset)) == 1)
	stopifnot(nrow(eset) == 92)
	
	stopifnot(sum(grepl(";",fData(eset)$proteinName)) == 2)
	stopifnot("sp|Q9Y6H1|CHCH2_HUMAN;sp|Q5T1J5|CHCH9_HUMAN" %in% fData(eset)$proteinName)
	
	cat(" --- testParseProgenesisPeptideMeasurementCsv: PASS ALL TEST  --- \n")
	
}


testParseMaxQuantProteinGroupTxt <- function(){
	
	cat(" --- testParseMaxQuantProteinGroupTxt:  --- \n")
	
	expDesign <- data.frame(row.names=1:20,isControl=c(rep(T,5),rep(F,15)),condition=sort(rep(paste("condition",1:4,sep="_"),5)))
	eset <- parseMaxQuantProteinGroupTxt(maxQuantProteinFileTxt,expDesign=expDesign, method="auc")
	
	stopifnot(ncol(fData(eset)) == 7)
	stopifnot(ncol(eset) == 20)
	stopifnot( nrow(parseMaxQuantProteinGroupTxt(maxQuantProteinFileTxt,expDesign=expDesign, method="spc")) == 996 )

	cat(" --- testParseMaxQuantProteinGroupTxt: PASS ALL TEST  --- \n")
	
}

testParseScaffoldPTMReport <- function(){
	
	cat(" --- testParseScaffoldPTMReport:  --- \n")
	df <- parseScaffoldPTMReport(scaffoldPtmReportFile1)
	stopifnot(nrow(df) == 883 )
	cat(" --- testParseScaffoldPTMReport: PASS ALL TEST  --- \n")

}


testAddScaffoldPTMFAnnotations <- function(){
	
	expDesignTMTTenPlex <- data.frame(condition=sort(rep(c(1:5),2)),isControl=c(T,T,rep(F,8)) )
	eset <- parseScaffoldRawFile(scaffoldPtmTMTRawDataFile1,expDesign = expDesignTMTTenPlex)
	
	cat(" --- testAddScaffoldPTMFAnnotations:  --- \n")
	eset <- addScaffoldPTMFAnnotations(eset,scaffoldPtmReportFile1)
	
	stopifnot(all( c("ptmPeptide","ptm","ptmLocProb","idScore","ptmLocMascotConfidence","pMassError") %in% names(fData(eset))))
	
	cat(" --- testAddScaffoldPTMFAnnotations: PASS ALL TEST  --- \n")
	
}

### TEST FUNCTIONS END


### TESTS
testGetFileType()

# progenesis
testGetProgenesisCsvExpressionColIndices()
testGetExpDesignProgenesisCsv()
testParseProgenesisProteinCsv()
testParseProgenesisFeatureCsv()
testParseProgenesisPeptideMeasurementCsv()

# scaffold 
testGetSkipLineNb()
testParseScaffoldRawFile()
testGetNbPlex()
testParseScaffoldPTMReport()
testAddScaffoldPTMFAnnotations()

# max quant
testParseMaxQuantProteinGroupTxt()

### TESTS END

