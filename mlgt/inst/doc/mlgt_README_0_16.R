### R code from vignette source '../mlgt_README_0_16'

###################################################
### code chunk number 1: Install seqinr (eval = FALSE)
###################################################
## install.packages("seqinr")


###################################################
### code chunk number 2: Install mlgt (eval = FALSE)
###################################################
## install.packages("mlgt_0.15.zip", repos=NULL)


###################################################
### code chunk number 3: Load library
###################################################
library(mlgt)


###################################################
### code chunk number 4: Set paths to external apps
###################################################
Sys.setenv(BLASTALL_PATH="C:/Users/Public/Apps/Blast/bin/blastall.exe",
		FORMATDB_PATH="C:/Users/Public/Apps/Blast/bin/formatdb.exe",
		FASTACMD_PATH="C:/Users/Public/Apps/Blast/bin/fastacmd.exe",
		MUSCLE_PATH="C:/Users/Public/Apps/Muscle/muscle3.8.31_i86win32.exe")


###################################################
### code chunk number 5: Set working directory (eval = FALSE)
###################################################
## analysisDir <-  "C:/Users/me/genoProject1/run1/analysis/"
## setwd(analysisDir)


###################################################
### code chunk number 6: Describe run
###################################################
system.file("namedBarcodes.fasta", package="mlgt")

# Load MIDs used to mark samples
fTagList <- read.fasta(system.file("namedBarcodes.fasta", package="mlgt"), 
			as.string=T) 
# Optionally, rename the barcodes to the samples used in this run
sampleBarcodeTable <- read.delim(system.file("tableOfSampleBarcodeMapping.tab", 
		package="mlgt"), header=T)
names(fTagList) <- sampleBarcodeTable$sample[
			match(names(fTagList), sampleBarcodeTable$barcode)]
# here we're using the same tags at both ends of the amplicons.
rTagList <- fTagList
#The names of the samples
sampleList <- names(fTagList)
# Load the marker sequences. 
myMarkerList <- read.fasta(system.file("HLA_namedMarkers.fasta", package="mlgt"),
			as.string=T)	
		
# The fasta file of sequence reads
inputDataFile <- system.file("sampleSequences.fasta", package="mlgt")


###################################################
### code chunk number 7: Prepare run
###################################################
# Creates object to store run settings
my.mlgt.Design <- prepareMlgtRun(projectName="myProject", 
				runName="myRun", samples=sampleList, 
				markers=myMarkerList, fTags=fTagList, 
				rTags=rTagList, inputFastaFile=inputDataFile, 
				overwrite="yes")


###################################################
### code chunk number 8: Inspect BLAST results
###################################################
# inspect BLAST results for a specific marker
thisMarker <- "DPA1_E2"
topHits <- getTopBlastHits(my.mlgt.Design@markerBlastResults)
#inspectBlastResults(topHits, thisMarker)


###################################################
### code chunk number 9: plot_BLAST_results
###################################################
inspectBlastResults(topHits, thisMarker)


###################################################
### code chunk number 10: print BLAST graphs to file (eval = FALSE)
###################################################
## # automatic output to pdf of blast result graphs for a list of markers.
## printBlastResultGraphs(my.mlgt.Design)


###################################################
### code chunk number 11: Run mlgt
###################################################
my.mlgt.Result <- mlgt(my.mlgt.Design)
save(my.mlgt.Result, file="thisRun.mlgtResult.Rdata")


###################################################
### code chunk number 12: Inspect run summary table
###################################################
my.mlgt.Result@runSummaryTable


###################################################
### code chunk number 13: Inspect results for a marker
###################################################
thisMarker <- "DPA1_E2"
my.mlgt.Result@markerSampleList[[thisMarker]]


###################################################
### code chunk number 14: Call genotypes
###################################################
my.genotypes <- callGenotypes(my.mlgt.Result)


###################################################
### code chunk number 15: Write genotype results table to file (eval = FALSE)
###################################################
## writeGenotypeCallsToFile(my.genotypes)


###################################################
### code chunk number 16: plot_genotype_evidence
###################################################
plotGenotypeEvidence(genotypeCall=my.genotypes[["DPA1_E2"]])	


###################################################
### code chunk number 17: Create a variant map from known alleles (eval = FALSE)
###################################################
## 
## markerImgtFileTable <- read.delim(system.file("marker.imgt.msf.list.tab", package="mlgt"),
##  				header=T)
## alignFilesSource <- 'ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/'
## # select a folder to store the alignments in. Here using current working directory.
## alignFilesDir <- getwd()	
## ## Download the allele alignments and create a 'variantMap' object for each marker and store them all in a list.
## knownAlleleDb <- list()
## for(thisMarker in names(myMarkerList)) {	
## 	fileName <-  markerImgtFileTable$imgtAlignFile[match(thisMarker, markerImgtFileTable$marker)]
## 	alleleAlignUrl  <- paste(alignFilesSource , fileName , sep="/")
## 	alleleAlignFile <- paste(alignFilesDir , fileName , sep="/")
## 	download.file(alleleAlignUrl,alleleAlignFile)
## 	knownAlleleDb[[thisMarker]] <- createKnownAlleleList(thisMarker,
## 		myMarkerList[[thisMarker]][1], alleleAlignFile)
## }


###################################################
### code chunk number 18: Call genotypes using known alleles (eval = FALSE)
###################################################
## my.genotypes <- callGenotypes(my.mlgt.Result,  mapAlleles=TRUE,
## 		alleleDb=knownAlleleDb)


###################################################
### code chunk number 19: Align_Report_profile
###################################################
alignReport(my.mlgt.Result,markers="DPA1_E2", samples="Sample-8", method="profile")


###################################################
### code chunk number 20: Align_Report_hist
###################################################
alignReport(my.mlgt.Result,markers="DPA1_E2", samples="Sample-8", method="hist")


###################################################
### code chunk number 21: Error_Correct
###################################################
my.mlgt.Result.Corrected <- errorCorrect(my.mlgt.Result)
# Produce an alignment report for the un-corrected and corrected results.
alignReport(my.mlgt.Result, method="profile", fileName="alignReport_my.mlgt.Result")
alignReport(my.mlgt.Result.Corrected, method="profile", fileName="alignReport_my.mlgt.Result.Corrected")


###################################################
### code chunk number 22: Combine_Results
###################################################
my.design.list <- list()
my.design.list[['A']] <- my.mlgt.Design
my.design.list[['A']]@samples <- sampleList[1:5]
my.design.list[['B']] <- my.mlgt.Design
my.design.list[['B']]@samples <- sampleList[6:10]
my.result.list <- lapply(my.design.list, FUN=function(x) mlgt(x))
my.result.list
combined.result <- combineMlgtResults(my.result.list)
combined.result


###################################################
### code chunk number 23: Show_list_use
###################################################
# Create a list of mlgtDesign objects, each with only one marker.
my.design.list <- list()
for(thisMarker in names(myMarkerList))  {
	my.design.list[[thisMarker]] <- my.mlgt.Design
	my.design.list[[thisMarker]]@markers <- myMarkerList[thisMarker]
	
}
# Use lapply to run mlgt() on each member of the list. 
# N.B. we are using errorCorrection within mlgt(), which slows it down a bit.
system.time(
	my.result.list <- lapply(my.design.list, 
			FUN=function(x) mlgt(x, errorCorrect=TRUE))
)


###################################################
### code chunk number 24: Prep_snowfall
###################################################
#install.packages('snowfall')
library(snowfall)

sfInit(parallel=TRUE, cpus=4, type="SOCK")	# set your number of processors here.
sfExport(list=ls())	# is this necessary?
sfLibrary(mlgt)		# the 'nodes' need to load a copy of the relevant libraries
sfLibrary(seqinr)	# is this one necessary?
# Then we run mlgt over the list of mlgtDesign objects. 
# Note that extra parameters can be passed to sfLapply().
system.time(
	sf.result.list <- sfLapply(my.design.list, mlgt, errorCorrect=TRUE)
)


###################################################
### code chunk number 25: Run_snowfall
###################################################
project.mlgt.Results <- combineMlgtResults(sf.result.list)


###################################################
### code chunk number 26: Custom_Call
###################################################
callGenotypes.custom <- function(table, maxPropUniqueVars=0.5) {
	table$status <- "notCalled"
	table$propUniqueVars <- table$numbVar/table$numbSeq
	table$status <- ifelse(table$propUniqueVars <= maxPropUniqueVars,"good", "bad")
	return(table)
}
my.custom.Genotypes <- callGenotypes(my.mlgt.Result, method="callGenotypes.custom")


###################################################
### code chunk number 27: Output sequences (eval = FALSE)
###################################################
## dumpVariantMap.mlgtResult(my.mlgt.Result)


###################################################
### code chunk number 28: Output all sequences (eval = FALSE)
###################################################
## dumpVariants(my.mlgt.Result)


###################################################
### code chunk number 29: Session info
###################################################
sessionInfo()


