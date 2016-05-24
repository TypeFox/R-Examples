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

proteinSeq1 <- "MSAGSSCSQTPSRAIPTRRVALGDGVQLPPGDYSTTPGGTLFSTTPGGTRIIYDRKFLMECRNSPVAKTPPKDLPAIPGVTSPTSDEPPMQASQSQLPSSPEDKRAGGEESQFEMDI"
proteinSeq2 <- "MVKKSRRRGAAQWAAVRAQAGLTATDENEDDLGLPPSPGDSSYYQDQVDEFHEARSRAVLAKGWNEVESGEEDGDEEEE"
proteinSeq3 <- "MSERMRPVVVDLPTSASSSMKVNG"	
proteinSeq4 <- "RKR"

### read protein db
proteinDB <- read.fasta(fastaFile,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

### INIT END

### TEST FUNCTIONS

testIsCon <- function(){
	cat("--- testIsCon: --- \n")
	stopifnot(sum(isCon(fData(eset)$proteinName)) == 0)
	cat("--- testIsCon: PASS ALL TEST --- \n")
}

testIsDecoy <- function(){
	cat("--- testIsDecoy: --- \n")
	stopifnot(sum(isDecoy(fData(eset)$proteinName)) == 198)
	cat("--- testIsDecoy: PASS ALL TEST --- \n")
}

testGetIdQvals <- function(){
	cat("--- testGetIdQvals: --- \n")
	stopifnot(max(getIdLevelQvals(fData(eset)$idScore, isDecoy(fData(eset)$proteinName))) < 1)
	cat("--- testGetIdQvals: PASS ALL TEST --- \n")
}


testAddIdQvalues <- function(){
	
	cat("--- testAddIdQvalues: --- \n")	
	stopifnot(all.equal(round(max(fData(addIdQvalues(eset))$idQValue ),2),0.28))
	stopifnot(all.equal(round(max(fData(addIdQvalues(rollUp(eset, featureDataColumnName = "proteinName"  )))$idQValue),2),0.29))
	cat("--- testAddIdQvalues: PASS ALL TEST --- \n")
	
	#progenesisPeptideCsvFile3 <- "/Users/ahrnee-adm/dev/R/workspace/SafeQuant/inst/testData/PeptidesSQAnalysis/peptides5.csv"
	#d <- parseProgenesisPeptideCsv(file=progenesisPeptideCsvFile3,expDesign=getExpDesignProgenesisCsv(progenesisPeptideCsvFile3))
	#e <- rollUp(d, method="sum", isProgressBar=T, featureDataColumnName= c("proteinName"))
	#x <- addIdQvalues(e)
	#
	#fData(x)$idQValue
	#isPtm <- nchar(as.character(fData(x)$ptm)) > 0
	#plot(sort(fData(x)$idScore[isPtm]),sort(fData(x)$idQValue[isPtm],decreasing=T), type="l" )
	#lines(sort(fData(x)$idScore[!isPtm]),sort(fData(x)$idQValue[!isPtm],decreasing=T), col=2)
	#
	#
	#d <- parseProgenesisProteinCsv(file=progenesisProteinCsvFile1,expDesign=getExpDesignProgenesisCsv(progenesisProteinCsvFile1))
	#x <- addIdQvalues(d)
	#fData(x)$idQValue
	#plot(sort(fData(x)$idScore),sort(fData(x)$idQValue,decreasing=T), type="l" )
	#
	#hist(fData(x)$idQValue)
	
	
}

testGetScoreCutOff <- function(){
	cat("--- testGetScoreCutOff: --- \n")
	stopifnot(all.equal(round(getScoreCutOff(fData(eset)$idScore, isDecoy(fData(eset)$proteinName)),2), 12.07))
	cat("--- testGetScoreCutOff: PASS ALL TEST --- \n")
}

testGetModifProteinCoordinates <- function(){
	
	ptm <- "[N-term] Acetyl (Protein N-term)|[6] Oxidation (M)"
	peptide <- "SSDAEMAVFGEAAPYLR"
	protein <-  "MPEPTIDESSDAEMAVFGEAAPYLRKSEKERIEAQNKPFDAK"
	
	#debug(getModifProteinCoordinates)
	
	cat("--- testGetModifProteinCoordinates: --- \n")
	mc <- getModifProteinCoordinates(ptm,peptide,protein)
	stopifnot(9 == mc[1])
	stopifnot(14 == mc[2])
	stopifnot(all.equal(c(34,49),getModifProteinCoordinates(fData(eset)$ptm[1],fData(eset)$peptide[1],proteinSeq1))) 
	stopifnot(37 == getModifProteinCoordinates(fData(eset)$ptm[2],fData(eset)$peptide[2],proteinSeq2)) 
	
	#### SCAFFOLD PTM FORMAT
	proteinSeq <-  "MMMMMMMMMMETPSPRPPPMRHRSSRSP"		
	peptideSeq <- "ETPSPRPPPMR"
	modifAnnot <- "T2 Phospho, S4 Phospho, M10 Oxidation"
	stopifnot(all.equal(c(12,14,20),  getModifProteinCoordinates(modifAnnot,peptideSeq,proteinSeq, format=2)))
	
	
	cat("--- testGetModifProteinCoordinates: PASS ALL TEST --- \n")
}

testGetMotifX <- function(){
	
	cat("--- testGetMotifX: --- \n")
	
	modifPos1 <- getModifProteinCoordinates(fData(eset)$ptm[1],fData(eset)$peptide[1],proteinSeq1)
	modifPos2 <- getModifProteinCoordinates(fData(eset)$ptm[2],fData(eset)$peptide[2],proteinSeq2)
	
	stopifnot("PGDYS*TTPG" == getMotifX(modifPos1,fData(eset)$peptide[1],proteinSeq1,4)[1])
	stopifnot("GLPPS*PGDS" == getMotifX(modifPos2,fData(eset)$peptide[2],proteinSeq2,4))
	
	#stopifnot("PGDYS*TTPG" == getMotifX(fData(eset)$ptm[1],fData(eset)$peptide[1],proteinSeq1,4)[1])
	#stopifnot("GLPPS*PGDS" == getMotifX(fData(eset)$ptm[2],fData(eset)$peptide[2],proteinSeq2,4))
	
	cat("--- testGetMotifX: PASS ALL TEST --- \n")
	
}

testAddPTMCoord <- function(){
	
	cat("--- testAddPTMCoord: --- \n")	
	eset <- .addPTMCoord(eset,proteinDB,motifLength = 4)
	stopifnot( "PGDYS*TTPG,TPGGT*RIIY"  ==  fData(eset)$motifX[1] )
	cat("--- testAddPTMCoord: PASS ALL TEST --- \n")
	
}


testSetFilter <- function(){
	
	### peptide data example
	
	idQValueThrs <- 0.01
	pMassTol <- c(-1,1)
	ptmRegExp <- "?"
	proteinNameRegExp <- "?"
	
	filter <- data.frame(
			
			fData(addIdQvalues(eset))$idQValue >= idQValueThrs # id score
			,(fData(eset)$pMassError < pMassTol[1]) | (fData(eset)$pMassError > pMassTol[2]) # precursor mass tolerance
			,isDecoy(fData(eset)$proteinName)	# decoy
			,isCon(fData(eset)$proteinName)	# contaminants	
			,!(regexpr(proteinNameRegExp,fData(eset)$proteinName,ignore.case=T) == 1) # protein ac 
			,!(regexpr(ptmRegExp,fData((eset))$ptm,ignore.case=T) == 1 ) # ptm 
	
	)
	
	cat("--- testSetFilter: --- \n")	
	
	stopifnot(sum(fData(.setFilter(eset,filter=filter))$isFiltered) == 230)
	stopifnot(sum(fData(.setFilter(eset,filter=cbind(filter,rep(T,nrow(eset)))))$isFiltered) == 900)
	
	cat("--- testSetFilter: PASS ALL TEST --- \n")
	
}

testGetPeptides <- function(){
	
	stopifnot(length(getPeptides(proteinSeq3)) == 3)
	stopifnot(paste(getPeptides(proteinSeq3),collapse="") == proteinSeq3)
	stopifnot(length(getPeptides(proteinSeq3,proteaseRegExp=.getProteaseRegExp("lys-c"))) == 2)

	stopifnot(length(getPeptides(proteinSeq4,nbMiscleavages=0)) == 3)
	stopifnot(length(getPeptides(proteinSeq4,nbMiscleavages=1)) == 5)
	stopifnot(length(getPeptides(proteinSeq4,nbMiscleavages=2)) == 6)
	
	### digest whole db
	if(F){

		digestedDB <- list()
		
		for(i in 1:length(names(proteinDB))){
		#for(i in 1:100){
			
			cat(i,"/",length(proteinDB),"\n")
			ac <- names(proteinDB)[i]
			
			for(peptide in getPeptides(proteinDB[[ac]])){
				
				digestedDB[[peptide]] <- c(digestedDB[[peptide]],ac)
				
			}
			
		}
	}
}

testGetNbDetectablePeptides <- function(){
	
	cat("--- testGetNbDetectablePeptides: --- \n")	
	stopifnot(getNbDetectablePeptides(getPeptides(proteinSeq1)) == 8)
	stopifnot(getNbDetectablePeptides(getPeptides(proteinSeq2),peptideLength=c(-Inf,Inf)) == length(getPeptides(proteinSeq2)))
	cat("--- testGetNbDetectablePeptides: PASS ALL TEST --- \n")
	
}

testGetNbMisCleavages <- function(){
	
	cat("--- testGetNbMisCleavages: --- \n")	
	peptide <- c("PEPTIDEK","PERTIDEK","PERKTIDEK","PERTIDE","RRPERPKK")
	protease <- "trypsin"
	stopifnot( all.equal(c(0,1,2,1,2),getNbMisCleavages(peptide)))
	cat("--- testGetNbMisCleavages:  PASS ALL TEST --- \n")	
	
}

testGetPeptidePerProtein <- function(){
	
	cat("--- testGetPeptidePerProtein: --- \n")
	stopifnot(getNbPeptidesPerProtein(eset)["prot52"] == 3)
	cat("--- testGetPeptidePerProtein: PASS ALL TEST --- \n")
	
}

testSetPeptidePerProtein <- function(){
	
	cat("--- testSetPeptidePerProtein: --- \n")
	stopifnot("nbPeptides" %in% names(fData(setNbPeptidesPerProtein(eset))))
	cat("--- testSetPeptidePerProtein: PASS ALL TEST --- \n")
	
}

testGetMeanCenteredRange <- function(){
	
	cat("--- testGetMeanCenteredRange: --- \n")
	stopifnot( round(mean(getMeanCenteredRange(fData(eset)$pMassError)),3) == round(mean(fData(eset)$pMassError),3) )
	cat("--- testGetMeanCenteredRange: PASS ALL TEST --- \n")
	
}

testIsStrippedACs <- function(){
	
	cat("--- testIsStrippedACs: --- \n")
	stopifnot(!isStrippedACs( sample(names(proteinDB),100) ) )
	stopifnot(isStrippedACs(fData(eset)$proteinName))
	cat("--- testIsStrippedACs: PASS ALL TEST --- \n")
}

testStripACs <- function(){
	
	cat("--- testStripACs: --- \n")
	stopifnot(isStrippedACs(stripACs(sample(names(proteinDB),100))))
	cat("--- testStripACs: PASS ALL TEST --- \n")
}

### TEST FUNCTIONS END

testGetAAProteinCoordinates <- function(){
	
	cat("--- testGetAAProteinCoordinates: --- \n")
	peptide <- "SSDAEMAVFGEAAPYLR"
	protein <-  "MPEPTIDESSDAEMAVFGEAAPYLRKSEKERIEAQNKPFDAK"
	stopifnot(length(getAAProteinCoordinates("SSDAEMAVFGEAAPYLR","MPEPTIDESSDAEMAVFGEAAPYLRKSEKERIEAQNKPFDAK","S")) == 2)
	stopifnot(getAAProteinCoordinates("SSDAEMAVFGEAAPYLR","MPEPTIDESSDAEMAVFGEAAPYLRKSEKERIEAQNKPFDAK","Y") == 23)
	cat("--- testGetAAProteinCoordinates: PASS ALL TEST --- \n")	
	
}

### TESTS

testIsCon()
testIsDecoy()
testGetIdQvals()
testAddIdQvalues()
testGetScoreCutOff()
testGetModifProteinCoordinates()
testGetMotifX()
testAddPTMCoord()
testSetFilter()
testGetPeptides()
testGetNbDetectablePeptides()
testGetNbMisCleavages()
testGetPeptidePerProtein()
testSetPeptidePerProtein()
testGetMeanCenteredRange()
testIsStrippedACs()
testStripACs()
testGetAAProteinCoordinates()

### TESTS END




