# TODO: Add comment
# 
# Author: ahrnee-adm
###############################################################################

### load / source
library("limma")
library(gplots) # volcano plot
library(seqinr)
library(optparse)
library(data.table)
library(epiR)
library(corrplot)
library(Biobase)

### INIT
if(!grepl("SafeQuant\\.Rcheck",getwd())){ # DEV mode 
	wd <- dirname(sys.frame(1)$ofile)
	setwd(dirname(sys.frame(1)$ofile))
	sqRootDir <- dirname(getwd())
	
	source(paste(sqRootDir,"/R/ExpressionAnalysis.R",collapse="",sep=""))
	source(paste(sqRootDir,"/R/SafeQuantAnalysis.R",collapse="",sep=""))
	source(paste(sqRootDir,"/R/Graphics.R",collapse="",sep=""))
	source(paste(sqRootDir,"/R/IdentificationAnalysis.R",collapse="",sep=""))
	source(paste(sqRootDir,"/R/Parser.R",collapse="",sep=""))
	source(paste(sqRootDir,"/R/TMT.R",collapse="",sep=""))
	source(paste(sqRootDir,"/R/UserOptions.R",collapse="",sep=""))
	
	source(paste(sqRootDir,"/R/Targeted.R",collapse="",sep=""))
	
	
}else{ # CHECK mode
	### wd already set to tests when running CHECK
	library(SafeQuant)
}



### INIT
### VARIOUS TEST FILES

# progenesis
progenesisFeatureCsvFile1 <- "testData/progenesis_feature_export1.csv"
progenesisPeptideMeasurementCsvFile1 <- "testData/progenesis_pep_measurement1.csv"

#progenesisPeptideMeasurementCsvFile1 <- "testData/tmp.csv"

progenesisProteinCsvFile1 <- "testData//progenesis_protein_export1.csv"
progenesisPeptideMeasurementFractionatedCsvFile1 <- "testData/progenesis_pep_measurement_fractionated1.csv"

progenesisProteinCsvFile2 <- "testData/2014/proteins2.csv"
progenesisFeatureCsvFile2 <- "testData/2014/peptides2.csv"


# scaffold
scaffoldTmt6PlexRawTestFile <- "testData/scaffold_tmt6plex_raw.xls"
scaffoldTmt10PlexRawTestFile <- "testData/scaffold_tmt10plex_raw.xls"

#scaffoldPtmTMTRawDataFile1 <- "testData/scaffoldPTM/Christoph-LE-Human-pH10fraction-TMT-20150630/Raw Data Report for Christoph-LE-Human-pH10fraction-TMT-20150630.xls"
#scaffoldPtmReportFile1 <- "testData/scaffoldPTM/Christoph-LE-Human-pH10fraction-TMT-20150630/Spectrum Report of Scaffold_PTM_P-TMT-pH10 Experiment.xls"
scaffoldPtmTMTRawDataFile1 <- "testData/scaffold_tmt10plex_raw_phospho.xls"
scaffoldPtmReportFile1 <- "testData/scaffoldPtm_spectrum_report.xls"

# maxquant
maxQuantProteinFileTxt <- "testData/maxquant_protein_groups.csv"

# db
fastaFile <- "testData/mouse_proteins.fasta"

### INIT END

## CREATE TEST DATA

set.seed(1234)
nbFeatures <- 900

peptide <- paste("pep",1:nbFeatures,sep="")
peptide[1] <- "VALGDGVQLPPGDYSTTPGGTLFSTTPGGTR"
peptide[2] <- "AQAGLTATDENEDDLGLPPSPGDSSYYQDQVDEFHEAR"

proteinName <- sort(rep(paste("prot",1:(nbFeatures/3),sep=""),3))
proteinName[1:200] <- paste("REV_",proteinName[1:200] ,sep="")
proteinName[1] <- "sp|Q60876|4EBP1_MOUSE"
proteinName[2] <- "sp|Q9JI13|SAS10_MOUSE"

idScore <- rep(0,length(proteinName))
idScore[1:200] <-rnorm(200,10,1)
idScore[c(1:2,201:900)] <- rnorm(702,15,1)

ptm <- rep("",900)
ptm[1] <- "[15] Phospho (ST)|[30] Phospho (ST)"
ptm[2] <- "[20] Phospho (ST)"

pMassError <- c(rnorm(200,0,1.5),rnorm(700,0,0.5))

charge <- round(runif(length(ptm),1.5,3.7))

peptideName <- paste(peptide,ptm)

proteinDescription <- sort(rep(paste("protDescription",1:(nbFeatures/3),sep=""),3))
isNormAnchor <- rep(T,nbFeatures)
isFiltered <- rep(F,nbFeatures)

m <- as.matrix( data.frame(rnorm(nbFeatures,1001),rnorm(nbFeatures,1001),rnorm(nbFeatures,1002),rnorm(nbFeatures,1002),rnorm(nbFeatures,1000),rnorm(nbFeatures,1000)) )
rownames(m) <- peptideName
colnames(m) <- c("A_rep_1","A_rep_2","B_rep_1","B_rep_2","C_rep_1","C_rep_2")

### phenoData: stores expDesign
#condition isControl
#A_rep_1         A     FALSE
#A_rep_2         A     FALSE
#B_rep_1         B     FALSE
#B_rep_2         B     FALSE
#C_rep_1         C      TRUE
#C_rep_2         C      TRUE

expDesign <- data.frame(condition=c("A","A","B","B","C","C"),isControl=c(F,F,F,F,T,T),row.names=colnames(m))
#expDesign <- data.frame(condition=c("A","A","B","B","C","C"),row.names=colnames(m))

featureAnnotations <- data.frame(
		 peptide
		, charge		 	
		,proteinName
		,proteinDescription
		,idScore
		,ptm
		,pMassError
		,isNormAnchor
		,isFiltered
		,row.names=peptideName)

eset <- createExpressionDataset(expressionMatrix=m,expDesign=expDesign,featureAnnotations=featureAnnotations)
sqa <- safeQuantAnalysis(eset)

# ABS. QUANT SIM. DATA 
cpc <- rep(2^(1:5),10)
set.seed(1234)
signal <- rnorm(length(cpc),cpc,cpc/10)
absEstSimData <- data.frame(cpc =  log10(cpc),signal = log10(signal))
absEstSimDataFit <- lm(cpc ~ signal, data=absEstSimData )


#data(proteomeMixLFQ,package="SafeQuant")
#data(proteomeMixTMT6,package="SafeQuant")
