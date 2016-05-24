# TODO: Add comment
# 
# Author: erikahrne
###############################################################################



### CMD OPTIONS
#' Command Line Option List
#' @export
option_list <- list(
		
		
### I/O
		make_option(c("-i", "--inputFile"), type="character", default="",
				help="I/O:  Input file: Progenesis (Feature,Protein or Peptide) .csv,
			or Scaffold Q+ (Raw Export, for TMT quant) .xls (REQUIRED)",
		),	
		make_option(c("-o", "--outputDir"), type="character", default=NA,
				help="I/O:  Results Output Directory [default FOLDER OF INPUTFILE]",
		),
		
		make_option(c("-l", "--resultsFileLabel"), type="character", default="SQ_Results",
				help="I/O: results file directory [default %default]", 
		),

		make_option(c("-f", "--fastaFile"), type="character", default="",
				help="I/O:  Protein DB .fasta file [default ./]",
		),
		
		make_option(c("-p", "--scaffoldPTMSpectrumReportFile"), type="character", default="",
				help="I/O:  Scaffold PTM Spectrum Report File [default ./]",
		),
		
### I/O END
		
# FILTER (--F)
		make_option(c("--FProteinAccessionSelection"), type="character", default=".",
				help="FILTER: --FP Filter features by Accession Regular Expression [default %default] (all features kept)",
				metavar="Protein Accession Reg. expr."),
		
		#### peptide analysis specfic
		make_option(c("--FModificationSelection"), type="character", default="",
				help="FILTER (LFQ PEP ONLY): --FM Only keep Peptides with modifications matching Regular Expression [default %default]
				 (all features kept). Peptide analysis ONLY",
				metavar="modification name Reg. expr."),
		
		make_option(c("--FFdrCutoff"), type="double", default=0.01,
				help="FILTER (LFQ ONLY): --FF Identification level False Discovery Rate Cutoff.  [0-1] [default %default]",
				metavar="Peptide/Protein FDR cutoff"),
		
		make_option(c("--FCoefficientOfVarianceMax"), type="double", default=Inf,
				help="FILTER: --FC Do not include features with C.V. above this threshold in statistical 
				test for differential expression [default %default]",
				metavar="Coefficent of Variance cutoff"),
		
		#### peptide analysis specfic
		make_option(c("--FDeltaMassTolerancePrecursor"), type="character", default="AUTO SET",
				help="FILTER (LFQ PEP ONLY): --FD Precursor mass Error Range filter (ppm) [default %default].
				Peptide imports ONLY",
				metavar="Mass Range [x,y]"),
		
		#### protein analysis specfic
		make_option(c("--FNumberOfPeptidesPerProteinMin"), type="integer", default=1,
				help="FILTER: --FN Only include those proteins with at least x identified peptides [default %default]
				Protein analysis ONLY.",
				metavar="Number of peptides"),
		
		#### peptide analysis specfic
		make_option(c("--FSitesPerPeptide"), type="integer", default=99999,
				help="FILTER: --FS Max Nb. Modifications Per Peptide [default Inf]
						Peptide analysis ONLY.",
				metavar="Max Number of PTM sites Per Petptide"),
		
		#### peptide analysis specfic
		make_option(c("--FLengthPeptide"), type="integer", default=1,
				help="FILTER: --FL Min Peptide Length (Nb. AA's) [default Inf]
						Peptide analysis ONLY.",
				metavar="Min Peptide Length (>=)"),
		
		
# FILTER (--F) END	
		
# STATISTICS (--S)
	
	make_option(c("--SAnchorProtein"), type="character", default=".",
			help="STATISTICS: --SA Normalize Intensities by selected protein(s) Regular Expression
			 [default %default] (use all proteins).",
			metavar="Protein Accession Reg. expr."),
	
	make_option(c("--SPvalueInclude"), action="store_true", default=FALSE,
			help="STATISTICS: --SP output eBayes moderated t-statistic p-values [default %default]"),
	
	make_option(c("--SNonPairWiseStatTest"), action="store_true", default=FALSE,
			help="STATISTICS: --SN non pairwise eBayes moderated t-statistic p-values [default %default]"),
	
# STATISTICS (--S) END

# EXPERIMENTAL DESIGN (--E)

	make_option(c("--EXperimentalDesign"), type="character", default=NA,
			help='EXPERIMENTAL DESIGN: --EX "," seperated samples, ":" separated conditions 
					Example: 1,2,3:4,5,6 
					   condition1 (REF) : channel 1,2,3
					   condition2: channel 4,5,6
					Note: for 10-plex default is "1,4,7,10:2,5,8:3,6,9"
					[default %default]'), 
	
	make_option(c("--EProteinQuantOff"), action="store_false", default=TRUE,
			help='EXPERIMENTAL DESIGN: --EP Disable Protein Level Quantification [default %default]'),
	
# EXPERIMENTAL DESIGN (--E) END

# PDF-REPORT (--P) 
	make_option(c("--PRatioCutOff"), type="double", default=1,
		help="PDF-REPORT: --PR Intensity ratio (fold change) cut-off used for graphics export. >1 [default %default]",
		metavar="Intensity ratio cutoff"),	

	make_option(c("--PQvaueCutOff"), type="double", default=0.01,
			help="PDF-REPORT: --PQ Qvalue cut-off used for graphics. 
			High-lighting features with a qval < specified value. [0-1] [default %default]",
			metavar="Differential expression qvalue cutOff"),	
	
#	make_option(c("--PSelectedGraphics"), type="character", default="",
#			help="PDF-REPORT: --PS Excluded Graphics: give letter for each plot to exclude ex: --PS iv 
#					(creates all plots but intensity density plots & volcano plot)
#					experimental design (e)
#					peptide feature score distrib related plots (f)
#					intensity distibution plots (i)
#					volcano plots (v)
#					hierarchical clustering plots (h)
#					differential expression fdr plot (d)	
#					[default (all plots) %default]"),		
# PDF-REPORT (--P) END

## TSV-REPORT (--T)
#	
#	make_option(c("--TFastaFile"), type="character", default="",
#			help="TSV-REPORT (LFQ PEP): -TF Protein Fasta File used to extract Modification Site Coordinates [default None]",
#			metavar=".fasta file path"),	
## TSV-REPORT (--T) END

# ADDITIONAL-REPORTS (--A)
	make_option(c("--ARDataFile"), action="store_true", default=FALSE,
		help="ADDITIONAL-REPORTS: --AR Save R objects in 'label'.RData file [default %default]"),

#	make_option(c("--AProtease"), type="character", default="KR",
#			help="ADDITIONAL-REPORTS: --TP protease [default (trypsin) %default]
#					1) trypsin
#					2) lys-c
#			
#					Option considered for iBAQ normalization", 
#	),

	make_option(c("--AIbaq"), action="store_true", default=FALSE,
			help="ADDITIONAL-REPORTS (LFQ PROT): --AI creates .tsv output file
					including protein iBAQ values. [default %default]"),
	
	make_option(c("--ATop3"), action="store_true", default=FALSE,
			help="ADDITIONAL-REPORTS (LFQ PEP): --AT creates .tsv output file
					including protein top3 values. [default %default]"),
	
# ADDITIONAL-REPORTS (--A) END

# TEST (peptide analysis specific)
	make_option(c("-t", "--test"), action="store_true", default=FALSE,
			help="TEST: test option, include first 2000 entries only [default %default]
			Peptide analysis ONLY."),
# TEST END
	make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
			help="Print extra output [default %default]")
	)

	
#' Read User Specified Command Line Options
#' @param version Safequant version number
#' @return user options list
#' @import  optparse
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getUserOptions <- function(version=version){
	
	epilogue <- "Examples:
	Progenesis LFQ Protein Quant:
	>Rscript safeQuant.R -i /path/to/peptide_measurment.csv

	Progenesis LFQ Protein Quant (QE):
	>Rscript safeQuant.R -i /path/to/peptide_measurment.csv --FL 7

	Progenesis LFQ Phospho Quant:
	>Rscript safeQuant.R -i /path/to/peptide_measurment.csv -f /path/to/proteins.fasta --FM phospho --FS 3 --EP

	Scaffold Q+ TMT Protein Quant:
	>Rscript safeQuant.R -i /path/to/Raw_Export.xls --EX 1,2,3,4:5,6,7:8,9,10
	
	Scaffold Q+ TMT  PEPTIDE PTM Quant (PHOSHO):
	>Rscript safeQuant.R -i /path/to/Raw_Export.xls -p /path/to/Spectrum_Export_Scaffold_PTM.xls --EX 1,2,3,4:5,6,7:8,9,10 --FM phospho --FS 3 -f /path/to/proteins.fasta
	"
	
	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	cmdOpt <- parse_args(OptionParser( prog=paste("SafeQuant",version), option_list=option_list, epilogue=epilogue))
	
	### CMD OPTIONS END						
	
	### SET USER OPTIONS
	userOptions <- list()

### VERBOSE
	#VERBOSE: verbose
	userOptions$verbose <- cmdOpt$verbose
### VERBOSE	END
	
# I/O
	#I/O: progenesisFilePath
	userOptions$inputFile <- cmdOpt$inputFile
	if( userOptions$inputFile == "" | !file.exists(userOptions$inputFile)){
		cat("ERROR. Please specify input file.",userOptions$inputFile, "Not found!","\n")
		q(status=-1)
	}
	
	#I/O: resultsFileLabel
	userOptions$resultsFileLabel <- cmdOpt$resultsFileLabel
	
	#I/O: outputDir
	userOptions$outputDir <- cmdOpt$outputDir
	if(is.na(userOptions$outputDir)){ # see default
		userOptions$outputDir <- dirname(userOptions$inputFile) 
	}
	if(!file.exists(userOptions$outputDir) & userOptions$outputDir != "" ){
		cat("ERROR. No such directory",userOptions$outputDir,"\n")
		q(status=-1)
	}else{
		userOptions$outputDir <- file.path(userOptions$outputDir, userOptions$resultsFileLabel)
	}
	
	#I/O: proteinFastaFile
	userOptions$proteinFastaFile <- NA
	if(nchar(cmdOpt$fastaFile) > 0 ){
		### check if file exists
		if(file.exists(cmdOpt$fastaFile)){
			userOptions$proteinFastaFile <- cmdOpt$fastaFile
		}else{
			cat("ERROR. File does not exist",cmdOpt$fastaFile,"\n")
			q(status=-1)
		}				
	}
	
	#I/O: scaffoldPTMSpectrumReportFile
	userOptions$scaffoldPTMSpectrumReportFile <- NA
	if(nchar(cmdOpt$scaffoldPTMSpectrumReportFile) > 0 ){
		### check if file exists
		if(file.exists(cmdOpt$scaffoldPTMSpectrumReportFile)){
			userOptions$scaffoldPTMSpectrumReportFile <- cmdOpt$scaffoldPTMSpectrumReportFile
		}else{
			cat("ERROR. File does not exist",cmdOpt$scaffoldPTMSpectrumReportFile,"\n")
			q(status=-1)
		}				
	}
	
	
# I/O END
	
# FILTER (--F)

	#FILTER: selectedProteinName
	userOptions$selectedProteinName <- cmdOpt$FProteinAccessionSelection
	
	#FILTER: selectedModifName
	userOptions$selectedModifName <- cmdOpt$FModificationSelection
	
	#FILTER: fdrCutoff
	userOptions$fdrCutoff <- cmdOpt$FFdrCutoff
	if(is.na(userOptions$fdrCutoff) | userOptions$fdrCutoff <= 0 | userOptions$fdrCutoff > 1 ){
		cat("ERROR. fdrCutOff must be in the range [0-1]. You specified",userOptions$fdrCutoff,"\n")
		q(status=-1)
	}
	
	#FILTER: precursorMassFilter
	if(cmdOpt$FDeltaMassTolerancePrecursor == "AUTO SET" ){
		userOptions$precursorMassFilter <- NA
	}else{
		### set by user -> add lower and upper mass bound to vector
		userOptions$precursorMassFilter <- gsub("(\\[)","",cmdOpt$FDeltaMassTolerancePrecursor)
		userOptions$precursorMassFilter <- gsub("(\\])","",userOptions$precursorMassFilter)
		userOptions$precursorMassFilter <- sort(as.numeric(unlist(strsplit(userOptions$precursorMassFilter,","))))
		
		### check input format precursorMassFilter 
		if((length(userOptions$precursorMassFilter) != 2) | sum(is.na(userOptions$precursorMassFilter)) > 0 ){
			cat("ERROR. Invalid FDeltaMassTolerancePrecursor", userOptions$minNbPeptidesPerProt, "\n")
			q(status=-1)
		}
	}
		
	#FILTER: cvCutOff
	userOptions$cvCutOff <- cmdOpt$FCoefficientOfVarianceMax
	if(is.na(userOptions$cvCutOff) | (userOptions$cvCutOff < 0)){
		print(paste("ERROR. cvCutOff must be > 0. You specified ", userOptions$cvCutOff))
		q(status=-1)
	}

	#FILTER: minNbPeptidesPerProt
	userOptions$minNbPeptidesPerProt <- cmdOpt$FNumberOfPeptidesPerProteinMin
	if(is.na(userOptions$minNbPeptidesPerProt) | userOptions$minNbPeptidesPerProt < 0 ){
		print(paste("ERROR. FNumberOfPeptidesPerProteinMin must be >= 0. You specified ", userOptions$minNbPeptidesPerProt))
		q(status=-1)
	}
	
	#FILTER: maxNbPTMsPerPeptide
	userOptions$maxNbPtmsPerPeptide <- cmdOpt$FSitesPerPeptide
	if(is.na(userOptions$maxNbPtmsPerPeptide) | userOptions$maxNbPtmsPerPeptide < 0 ){
		print(paste("ERROR. FSitesPerPeptide must be >= 0. You specified ", userOptions$maxNbPtmsPerPeptide))
		q(status=-1)
	}
	
	#FILTER: maxNbPTMsPerPeptide
	userOptions$minPeptideLength <- cmdOpt$FLengthPeptide
	if(is.na(userOptions$minPeptideLength) | userOptions$minPeptideLength < 0 ){
		print(paste("ERROR. FLengthPeptide must be >= 0. You specified ", userOptions$minPeptideLength))
		q(status=-1)
	}
	
# FILTER (--F) END

# STATISTICS
	
	#STATISTICS: normAC
	userOptions$normAC <- cmdOpt$SAnchorProtein

	#STATISTICS: eBayes
	userOptions$eBayes <- cmdOpt$SPvalueInclude
	
	#STATISTICS: SNonPairWiseStatTest
	userOptions$SNonPairWiseStatTest <- cmdOpt$SNonPairWiseStatTest
	
# STATISTICS END	
	
# EXPERIMENTAL DESIGN

	#EXPERIMENTAL DESIGN: EXperimentalDesign
	userOptions$expDesignTag <- cmdOpt$EXperimentalDesign
	
	userOptions$proteinQuant <- cmdOpt$EProteinQuant
	#userOptions$proteinQuant <- userOptions$selectedModifName != "."
		
# EXPERIMENTAL DESIGN END

# PDF-REPORT (--P)

	# PDF-REPORT: ratioCutOff
	userOptions$ratioCutOff <- cmdOpt$PRatioCutOff
	if(is.na(userOptions$ratioCutOff) | userOptions$ratioCutOff < 1){
		cat("ERROR. ratioCutoff must be > 1. You specified",userOptions$ratioCutOff,"\n")
		q(status=-1)
	}
	
	# PDF-REPORT: deFdrCutoff
	userOptions$deFdrCutoff <- cmdOpt$PQvaueCutOff
	if(is.na(userOptions$deFdrCutoff) | userOptions$deFdrCutoff <= 0 | userOptions$deFdrCutoff > 1 ){
		cat("ERROR. deFdrCutoff must be in the range [0-1]. You specified",userOptions$deFdrCutoff,"\n")
		q(status=-1)
	}

	# PDF-REPORT: PSelectedGraphics
#	userOptions$isDispExpDesign <- !regexpr("e",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isFdrPlots <- !regexpr("f",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isIntensityDistributionPlots <- !regexpr("i",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isVolcanoPlots <- !regexpr("v",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isHClustPlot <- !regexpr("h",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isDeFdrPlot <- !regexpr("d",cmdOpt$PSelectedGraphics) > -1


# PDF-REPORT (--P) END

# TSV-REPORT (--T)
#	
#	# TSV-REPORT: proteinFastaFile
#	userOptions$proteinFastaFile <- NA
#	if(nchar(cmdOpt$TFastaFile) > 0 ){
#		### check if file exists
#		if(file.exists(cmdOpt$TFastaFile)){
#			userOptions$proteinFastaFile <- cmdOpt$TFastaFile
#		}else{
#			cat("ERROR. File does not exist",cmdOpt$TFastaFile,"\n")
#			q(status=-1)
#		}				
#	}
#
#	# TSV-REPORT: proteaseTarget	(deprecated)
#	userOptions$protease <- cmdOpt$TProtease	

# TSV-REPORT (--T) END

# ADDITIONAL-REPORTS (--A)

	#ADDITIONAL-REPORTS iBaqTsvFile
	userOptions$iBAQ <- cmdOpt$AIbaq
#	userOptions$iBAQFile <- paste(userOptions$outputDir,userOptions$resultsFileLabel,"_iBAQ.tsv",sep="")

    #ADDITIONAL-REPORTS top3TsvFile
	userOptions$top3 <- cmdOpt$ATop3
#	userOptions$top3File <- paste(userOptions$outputDir,userOptions$resultsFileLabel,"_top3.tsv",sep="")

	#ADDITIONAL-REPORTS rDataFile, isSaveRObject	
	userOptions$isSaveRObject <- cmdOpt$ARDataFile
	#userOptions$rDataFile <- paste(userOptions$outputDir,userOptions$resultsFileLabel,".rData",sep="")
	
	
# ADDITIONAL-REPORTS (--A) END


# TEST

	### test run to define parameters (peptide analysis specific)
	userOptions$test <- cmdOpt$test

# TEST END
	return(userOptions)

}


#userInputTag <- "1,2,3:4,5,6"
# tag: 1,2:3:4,5,6 
#condition isControl
#1 Condition 1      TRUE
#2 Condition 1      TRUE
#3 Condition 1     TRUE
#4 Condition 2     FALSE
#5 Condition 2     FALSE
#6 Condition 2     FALSE

#' Create experimental design data.frame from user input string
#' @param tag tag
#' @param expDesignDefault data.frame  
#' @return data.frame describing experimental design
#' @export
#' @note  No note
#' @details  tag: 1,2:3:4,5,6 
#'		condition isControl
#'	1 Condition 1 TRUE
#'	2 Condition 1 TRUE
#'	3 Condition 1 TRUE
#'	4 Condition 2 FALSE
#'	5 Condition 2 FALSE
#'	6 Condition 2 FALSE
#' @references NA 
#' @examples print("No examples")
expDesignTagToExpDesign <- function(tag, expDesignDefault){
	
	sampleOrder <- as.numeric(unlist(strsplit(tag,"[\\,\\:]")))
	
	# make sure no duplicates, within range etc.
	if(is.na(sampleOrder[1]) 
			| (max(table(sampleOrder))>1) 
			| (min(sampleOrder) < 1)
			| max(as.numeric(sampleOrder)) > nrow(expDesignDefault)
			){
		stop("ERROR: getExpDesign, INVALID EXPERIMENTAL DESIGN ",tag,"\n")
		
	}
	
	expDesign <- data.frame(row.names=sampleOrder,condition=rep(NA,length(sampleOrder)), isControl=rep(FALSE,length(sampleOrder))  )
	
	condNb <- 1
	for(cond in unlist(strsplit(tag,":"))){
		
		#cat(as.character(unlist(strsplit(cond,","))), paste("Condition",condNb) , "\n")
		expDesign[as.character(unlist(strsplit(cond,","))),]$condition <- paste("Condition",condNb)
		condNb <- condNb + 1
	}       
	
	expDesign[ expDesign[,1] == "Condition 1" ,]$isControl <- T
	expDesign[,1] <- as.factor(expDesign[,1]) ### has to be factor and not character
	
	### get original sample names
	rownames(expDesign) <- rownames(expDesignDefault)[as.numeric(rownames(expDesign))]
	#expDesignUser$condition <- expDesign[rownames(expDesignUser) ,]$condition
	
	# get original condition names, unless conditions have been split
	
	# if more conditions than originally use cond_1, cond_n labellinf @TODO can be done better, to avoid loosing org condition names
	if(length(unique(expDesign$condition)) > length(unique(expDesignDefault$condition))) return(expDesign)
	
	# check if conditions are split in new expDesign
	for(cond in unique(expDesign$condition)){
		runs <- rownames(expDesign)[expDesign$condition == cond]
		
		if(length(unique(expDesignDefault[runs,]$condition)) > 1){
			return(expDesign)
			#stop("SPLIT")
		} 
	}
	
	# get original condition names
	expDesign$condition <- expDesignDefault[rownames(expDesign),]$condition
	
	return(expDesign)
	
}

