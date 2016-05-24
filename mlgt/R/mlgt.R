#' mlgt: Multi-locus geno-typing
#'
#' \tabular{ll}{
#' Package: \tab mlgt\cr
#' Type: \tab Package\cr
#' Version: \tab 0.16\cr
#' Date: \tab 2012-03-27\cr
#' Author: \tab Dave T. Gerrard <david.gerrard@@manchester.ac.uk>\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' mlgt sorts a batch of sequence by barcode and identity to 
#' templates. It makes use of external applications BLAST and 
#' MUSCLE. Genotypes are called and alleles can be compared
#' to a reference list of sequences. 
#' More information about each function can be found in 
#' its help documentation.
#'
#' Some text
#'
#' The main functions are: 
#' \code{\link{prepareMlgtRun}}, \code{\link{mlgt}},  
#' \code{\link{callGenotypes}}, \code{\link{createKnownAlleleList}},  
#' 
#' ...
#' 
#' @references BLAST - Altschul, S. F., W. Gish, W. Miller, E. W. Myers, and D. J. Lipman (1990). Basic local alignment search tool. Journal of molecular biology  215 (3), 403-410.
#' @references MUSCLE - Robert C. Edgar (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5), 1792-97.
#' @references IMGT/HLA database - Robinson J, Mistry K, McWilliam H, Lopez R, Parham P, Marsh SGE (2011) The IMGT/HLA Database. Nucleic Acids Research 39 Suppl 1:D1171-6

#' @import seqinr
#' @docType package
#' @name mlgt-package
NULL


# require(seqinr)

# n/b @ slot not currently recognised by roxygen.
# wanted reference to be of SeqFastadna, but unrecognised even with seqinr loaded.
#' An S4 class to hold all unique variants found/known for a marker.
#' 
#' @slot reference
#' @slot variantSource
#' @slot variantMap
#' @slot inputVariantCount
#' @slot uniqueSubVariantCount
setClass("variantMap", representation(reference='ANY', variantSource='character',variantMap='list',
							inputVariantCount='integer', uniqueSubVariantCount='integer'))


setClass("varCount", 
	representation(
		rawTotal="numeric",
		rawUniqueCount ="numeric",
		usedRawTotal="numeric",
		usedRawUniqueCount="numeric",
		subAlignTotal="numeric",
		subAlignUniqueCount="numeric", 
		varCountTable="data.frame"
	)
)


#' Create \code{\link{variantMap}} object from allele alignment
#' 
#' Create a \code{\link{variantMap}} object to store known alleles for a marker
#'
#' To compare variants produced using \code{\link{mlgt}} the sequences of the known
#' alleles must be aligned to the same marker used to find the variants. 
#' The resulting sub-sequence alignment may have identical sequences for different 
#' alleles. If that happens, those alleles are condensed into one and their names
#' concatenated.
#' User can supply files with marker sequences pre-aligned to the reference alleles. 
#' 
#' @param markerName A specific marker name
#' @param markerSeq something
#' @param alignedAlleleFile a sequence alignment 
#' @param alignFormat the format of alignedAlleleFile. "msf" (the default) or "fasta"
#' @param sourceName A character string to record the source of the alignment. Defaults to 
#'  the value of alignedAlleleFile
#' @param userAlignment The specified 'alignedAlleleFile' already includes the marker sequence. Default = FALSE.
#' 
#' @return a \code{\link{variantMap}} object named by markerName
#'
#' @export
#' @docType methods
createKnownAlleleList <- function(markerName, markerSeq, alignedAlleleFile, alignFormat="msf", sourceName=alignedAlleleFile, userAlignment=FALSE)  {
	## The aligned alleles must have unique names and the markerName must be different too. TODO: test for this.
	## TODO put default input file format (MSF). Make check for fasta format (can skip first part if already fasta).
	## clean up (remove) files. This function probably doesn't need to keep any files
	## Use defined class for return object giving marker sequence used as reference. 
	#alignedAlleles <- read.msf(alignedAlleleFile)
	
	musclePath <- Sys.getenv("MUSCLE_PATH")
	
	if(nchar(musclePath) < 1) {
			stop(paste("MUSCLE_PATH","has not been set!"))
	}

	switch(alignFormat,
		"msf" =   {
			alignedAlleles <- read.alignment(alignedAlleleFile, format="msf")
			alignedAlleleFastaFile <- paste(markerName, "alignedAlleles.fasta", sep=".")
			write.fasta(alignedAlleles[[3]], alignedAlleles[[2]], file.out=alignedAlleleFastaFile)
			}, 
		"fasta" = {
			alignedAlleleFastaFile = alignedAlleleFile
			},
		stop("Unrecognised alignment file format\n")
	)

	if(userAlignment) {	# user supplied alignment

		markerToAlleleDbAlign <- alignedAlleleFastaFile
	} else {	# align the marker to the aligned alleles.
			
		markerSeqFile <- paste(markerName, "markerSeq.fasta", sep=".")
		write.fasta(markerSeq, markerName, file.out=markerSeqFile )
		#muscle -profile -in1 existing_aln.afa -in2 new_seq.fa -out combined.afa
		markerToAlleleDbAlign <- paste(markerName, "allignedToAlleles.fasta", sep=".")
		# profile alignment of marker to existing allele alignment
		muscleCommand <- paste(musclePath, "-quiet -profile -in1", alignedAlleleFastaFile, "-in2", markerSeqFile, "-out" ,markerToAlleleDbAlign )
		system(muscleCommand)

	}

	### this section copied from getSubSeqsTable()  Could be recoded as function?

	# Extract portion corresponding to reference. 

	rawAlignment <- read.fasta(markerToAlleleDbAlign , as.string=T)		# do not use read.alignment() - broken
	alignedMarkerSeq <- s2c(rawAlignment[[markerName]])
	#alignedMarkerSeq <- s2c(rawAlignment[[thisMarker]])	# was this a bug?
	subStart <- min(grep("-",alignedMarkerSeq ,invert=T))
	subEnd <- max(grep("-",alignedMarkerSeq ,invert=T))
	alignedSubSeqs <- lapply(rawAlignment, FUN=function(x)	substr(x[1], subStart, subEnd))

	alignedSubTable <- data.frame(name=names(alignedSubSeqs ) , subSeq.aligned= as.character(unlist(alignedSubSeqs )))
	alignedSubTable$subSeq.stripped <-  gsub("-", "",alignedSubTable$subSeq.aligned )
	# remove marker sequence from allele subalignment list.
	alignedSubTable <- subset(alignedSubTable, name != markerName)
	# TODO: 2 issues here. 1. concatentated allele name lists get truncated at 500 chars (probably by newline insertion).
	# 2. blast is not returning the full name. 
	# HLA_A2 central region has 213 alleles all identical. The concatentated sequence IDs are 2030 chars long.
	# Could try using an '_' as this might prevent word boundary splitting in one or both cases.
	# Solves the first problem but still cannot blast with this. Need to also adjust blast usage.
	alleleMap <- split(as.character(alignedSubTable$name), alignedSubTable$subSeq.stripped)	# group allele nams sharing a common sequence.
	##alleleMap <- paste(unlist(split(as.character(alignedSubTable$name), alignedSubTable$subSeq.stripped)),collapse="|")
	# TODO: THis next line is what I want to do but it makes BLASTALL crash. Something to do with long sequence IDs in fasta files?
	alleleMap <- lapply(alleleMap, FUN=function(x) paste(x,collapse="_"))	# do not use '|' or '*' or ':'

	#if(remTempFiles) {
	#	file.remove(markerSeqFile)
	#}
	#return(list(reference=as.SeqFastadna(markerSeq, markerName), alleleMap=alleleMap, inputAlleleCount = length(unlist(alleleMap)), uniqueSubAlleleCount=length(alleleMap)))
	return(new("variantMap", reference=as.SeqFastadna(markerSeq, markerName), variantSource=sourceName, variantMap=alleleMap, inputVariantCount = length(unlist(alleleMap)), uniqueSubVariantCount=length(alleleMap)))
}




#' An S4 class that holds information about an mlgt analysis.
#'
#' Returned by \code{\link{prepareMlgtRun}}. Used as sole input for \code{\link{mlgt}}
#'
#' \describe{
#'   \item{projectName}{In which project does this run belong}
#'   \item{runName}{Which run was this. An identifier for the sequnece run}
#'   \item{markers}{A \emph{list} of named sequences.} 
#'   \item{samples}{A vector of sample names} 
#'   \item{fTags}{A vector of named sequence of MIDs used to barcode samples at the 5' end.}
#'   \item{rTags}{A vector of named sequence of MIDs used to barcode samples at the 3' end.}
#'   \item{inputFastaFile}{The name of the file containing sequences. Currently only fasta format is supported. It is up to you to pre-filter the sequences.}
#'  }
#'
#' @slot projectName In which project does this run belong
#' @slot runName Which run was this. An identifier for the sequnece run
#' @slot markers	A \emph{list} of named sequences. 
#' @slot samples A vector of sample names 
#' @slot fTags A vector of named sequence of MIDs used to barcode samples at the 5' end.
#' @slot rTags A vector of named sequence of MIDs used to barcode samples at the 3' end.
#' @slot inputFastaFile The name of the file containing sequences. Currently only fasta format is supported. It is up to you to pre-filter the sequences.
#' @seealso \code{\link{prepareMlgtRun}}, \code{\link{mlgt}}
setClass("mlgtDesign", 
	representation(
		projectName="character", 
		runName="character",
		markers="list",
		samples="character",
		fTags="list",
		rTags="list", 
		inputFastaFile="character", 
		markerBlastResults="character",
		fTagBlastResults="character",
		rTagBlastResults="character"
	)
)


setMethod("show", "mlgtDesign", definition= function(object="mlgtDesign"){
	cat("Design for mlgt run:\n")
	cat(paste("Project:",object@projectName,"\n"))
	cat(paste("Run:",object@runName,"\n"))
	cat(paste("Samples:",length(object@samples),"\n"))
	cat(paste("fTags:",length(object@fTags),"\n"))
	cat(paste("rTags:",length(object@rTags),"\n"))
	cat(paste("Markers:",length(object@markers),"\n"))
	}
)





#' An S4 class to hold results from \code{\link{mlgt}}
#'
#' Extends \code{\link{mlgtDesign}}
#'
#' \describe{
#'   \item{projectName}{In which project does this run belong}
#'   \item{runName}{Which run was this. An identifier for the sequnece run}
#'   \item{markers}{A \emph{list} of named sequences.} 
#'   \item{samples}{A vector of sample names} 
#'   \item{fTags}{A vector of named sequence of MIDs used to barcode samples at the 5' end.}
#'   \item{rTags}{A vector of named sequence of MIDs used to barcode samples at the 3' end. May be same as \code{fTags}}
#'   \item{inputFastaFile}{The name of the file containing sequences. Currently only fasta format is supported. It is up to you to pre-filter the sequences.}
#'   \item{runSummaryTable}{A summary table with one row per marker}
#'   \item{alleleDb}{A list of objects of class \code{\link{variantMap}}. Contains all variants returned by \code{\link{mlgt}}}
#'   \item{markerSampleList}{A list of tables, one table per marker giving results for each sample/MID}
#'  }
#'
#' @seealso \code{\link{mlgtDesign}}, \code{\link{prepareMlgtRun}}, \code{\link{mlgt}}
#'
setClass("mlgtResult", 
	representation(
			runSummaryTable="data.frame",
			alleleDb="list" ,
			markerSampleList="list",
			varCountTables="list"
	),
	contains="mlgtDesign"
)


setMethod("show", "mlgtResult", definition= function(object="mlgtResult"){
	cat("Results for mlgt run:\n")
	cat(paste("Project:",object@projectName,"\n"))
	cat(paste("Run:",object@runName,"\n"))
	cat(paste("Samples:",length(object@samples),"\n"))
	cat(paste("fTags:",length(object@fTags),"\n"))
	cat(paste("rTags:",length(object@rTags),"\n"))
	cat(paste("Markers:",length(object@markers),"\n"))
	print(object@runSummaryTable)
	}
)


#' Return top blast hits
#'
#' Auxillary function
#'
#' @param blastTableFile The name of a file of tabulated blast results.
#' @return A reduced blast table with one hit per query
#' @export
getTopBlastHits <- function(blastTableFile)  {		# returns the first hit for each query in the table. May now be partially redundant if selecting for number of blast hits returned..
	blastResults <- read.delim(blastTableFile, header=F)
	## Fields: 
	# Query id,Subject id,% identity,alignment length,mismatches,gap openings,q. start,q. end,s. start,s. end,e-value,bit score
	names(blastResults) <- c("query", "subject", "percent.id", "ali.length", "mismatches", "gap.openings", "q.start","q.end", "s.start","s.end", "e.value", "bit.score")
	topHits <- blastResults[match(unique(blastResults$query), blastResults$query),]
}


#INTERNAL. Need to document?
#' Align and trim sequences for marker/sample pair
#' 
#' A list of sequences mapped to both \option{thisMarker} and \option{thisSample} is created and these sequences are aligned to \option{markerSeq}.
#'
#' This internal function is called by \code{\link{mlgt}}
#' 
#' @param thisMarker A specific marker name
#' @param thisSample A specific sample name
#' @param sampleMap A list of sequence IDs assigned to each marker. Each element named by marker name.
#' @param fMarkerMap A list of sequence IDs assigned to each sample using BLAST hits in forward orientation. Each element named by sample name.
#' @param rMarkerMap A list of sequence IDs assigned to each sample using BLAST hits in reverse orientation. Each element named by sample name.
#' @param markerSeq The sequence of \option{thisMarker}
#' @param maxVarsToAlign If total assigned sequences exceeds 'minTotalCount', then only the 'maxVarsToAlign' most abundant variants are used.
#' @param minTotalCount How many assigned sequences to allow before limiting the number of raw variants to allign.
#' @param errorCorrect Use error correection on alignment of raw variants
#' @param correctThreshold Maximum proportion of raw reads at which (minor allele) bases and gaps are corrected.
#' @param minLength Reads below this length are excluded (they are very likely to be primer-dimers).
#'
#' @return A table of unique variants and their counts. The sequences have been trimmed to the portion aligned with \option{markerSeq}
#'
getSubSeqsTable <- function(thisMarker, thisSample, sampleMap, fMarkerMap,rMarkerMap, markerSeq,  maxVarsToAlign=30, minTotalCount=500, errorCorrect=FALSE, correctThreshold=0.01, minLength=70)  {
	if(exists("verbose")) cat("getSubSeqsTable",thisSample, thisMarker,"\n")
	# check paths to auxillary programs
	pathNames <- c("FASTACMD_PATH","MUSCLE_PATH")
	for(thisPath in pathNames)  {
		if(nchar(Sys.getenv(thisPath)) < 1) {
			stop(paste(thisPath,"has not been set!"))
		}
	}
	musclePath <- Sys.getenv("MUSCLE_PATH")
	fastacmdPath <- Sys.getenv("FASTACMD_PATH")
	if(length(grep(" ",musclePath, fixed=T))  > 0 ) musclePath <- shQuote(musclePath)
	if(length(grep(" ",fastacmdPath, fixed=T))  > 0 ) fastacmdPath <- shQuote(fastacmdPath)

	thisVarCount <- new("varCount", rawTotal=0,
			rawUniqueCount =0,
			usedRawTotal=0,
			usedRawUniqueCount=0,
			subAlignTotal=0,
			subAlignUniqueCount=0, 
			varCountTable=data.frame())
	varCountTable <- data.frame()
	#thisMarker <- "DPA1_E2"
	#thisSample <- "MID-1"
	fPairSeqList <- intersect(sampleMap[[thisSample]], fMarkerMap[[thisMarker]]) # fMarkerMap[[thisMarker]]
	rPairSeqList <- intersect(sampleMap[[thisSample]], rMarkerMap[[thisMarker]])

	# extract raw seqs from blastdb for forward hits. # THIS FUNCT CURRENTLY UNUSED BECAUSE NEED TO ASSESS IF ANYTHING TO QUERY BEFORE RUNNING fastacmd
	extractRawSeqsCommand <- function(idList,strand=1, fileName)  {
		queryList <- paste(idList , collapse=",")
		fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", fileName,  "-s", queryList)	
		return(fastacmdCommand)
	}

	fRawSeqs <- rRawSeqs <- list()
	
	## 12th Dec 11. Edited following rawseq extraction because failed when using large dataset. Replaced '-s' (queryList) option with '-i' (inputfile)
	if(length(fPairSeqList ) > 0)  {
		fIdFileName <- paste("test", thisMarker, thisSample, "fIdFile.txt",sep=".")
		write(fPairSeqList , file=fIdFileName )
		fRawSeqFileName <- paste("test", thisMarker, thisSample, "fRawSeqExtract.fasta",sep=".")
		strand <- 1
		fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", fRawSeqFileName,  "-i", fIdFileName)	
		system(fastacmdCommand)
		fRawSeqs <- read.fasta(fRawSeqFileName , as.string=T)
	}

	if(length(rPairSeqList ) > 0)  {
		rIdFileName <- paste("test", thisMarker, thisSample, "rIdFile.txt",sep=".")
		write(rPairSeqList , file=rIdFileName )		
		rRawSeqFileName <- paste("test", thisMarker, thisSample, "rRawSeqExtract.fasta",sep=".")
		strand <- 2
		fastacmdCommand  <- paste(fastacmdPath, "-p F -t T -d", "inputSeqs" , "-S", strand, "-o", rRawSeqFileName,  "-i", rIdFileName)
		system(fastacmdCommand)
		rRawSeqs <- read.fasta(rRawSeqFileName , as.string=T)
	}


	# Make file of unique seqs. Can name with sequence 

	### STRAND!!! - Done!

	# ver 0.14 March2012 - New algorithm required. Take minimum unique raw variants to achieve:-
	################ NO 50% (minPropReads) of reads 
	# or 30 most abundant unique variants (maxVarsToAlign)
	# or 500 reads (minTotalCount), whichever is larger. 

	rawSeqs <- c(fRawSeqs ,rRawSeqs )

	# filter out sequences shorter than minLength. Value of 0 means no filter.
	if(minLength > 0) rawSeqs <- rawSeqs[which(nchar(rawSeqs) > minLength)]


	totalRaw <- length(rawSeqs)

	if(totalRaw  < 1)  {
		#return(varCountTable)
		return(thisVarCount)
	}
	rawSeqCountTable <-  as.data.frame(table(unlist(rawSeqs)))
	rawVariantCount <- nrow(rawSeqCountTable)
	names(rawSeqCountTable) <- c("var", "count")
	rawSeqCountTable  <- rawSeqCountTable[order(rawSeqCountTable$count, decreasing=T),]	# most abundant raw vars at top of list.

	#minTotalCount <- ceiling(totalRaw * minPropReads)	# want to catch at least this many raw reads.
	#enoughReadsIndex <- cumsum(rawSeqCountTable$count)
	useIndex <- min(maxVarsToAlign, nrow(rawSeqCountTable))
	if(totalRaw > minTotalCount) {
		#use only most abundant 30 (maxVarsToAlign) unique variants.
		rawSeqCountTable <- rawSeqCountTable[1:useIndex,] 
	} 

	usedTotalRaw <- sum(rawSeqCountTable$count)
	usedVarCount <- nrow(rawSeqCountTable)
	cat(paste(thisSample,": Using", nrow(rawSeqCountTable), "variants, accounting for", usedTotalRaw, "of", totalRaw, "reads\n"))

	rawVariantFile <- paste("test", thisMarker, thisSample, "raw.variants.fasta",sep=".")  #"localVariants.fasta"
	#rawVariantFile <- paste(runPath, rawVariantFileName, sep="/")
	#v0.14# write.fasta(as.list(c(markerSeq,as.character(rawSeqCountTable$var))) ,c(thisMarker,as.character(rawSeqCountTable$var)),file.out=rawVariantFile ,open="w")	# align marker AND raw variants
	write.fasta(as.list(as.character(rawSeqCountTable$var)) ,as.character(rawSeqCountTable$var),file.out=rawVariantFile ,open="w") 	# align just raw variants


	# Align all seqs
	rawAlignFile <- paste("test", thisMarker, thisSample, "raw.align.fasta",sep=".")  #"localAlign.fasta"
	#rawAlignFile <- paste(runPath, rawAlignFileName, sep="/")
	##Removed for ver0.14 
	#if(rawVariantCount > 800) {
	#	muscleCommand <- paste(musclePath, "-in", rawVariantFile , "-out", rawAlignFile , "-diags -quiet -maxiters 2" )	# faster alignment
	#	warning(paste("Using fast MUSCLE alignment for ", thisMarker, thisSample, rawVariantCount, "sequences\n"))
	#} else {
		muscleCommand <- paste(musclePath, "-in", rawVariantFile , "-out", rawAlignFile , "-diags -quiet" )
	#}
	#cat(paste(muscleCommand, "\n"))
	system(muscleCommand)

	#v0.14# error Correct if required
	if(errorCorrect)  {		
			alignedVars <- read.fasta(rawAlignFile, as.string=T)
			seqCounts <- rawSeqCountTable$count[match(names(alignedVars),rawSeqCountTable$var)]
			#cat(seqCounts)
			#TODO MUST: strip off aligned marker and add back on - do not want this to be 'corrected' to match all other seqs. 
			#seqCounts <- varCountTable$count
			## important to 'unpack' the alignment so that each sequence occurs the correct number of times.
			sampleSeqs <- rep(alignedVars , seqCounts)
			thisAlign <- as.alignment(sum(seqCounts), names(sampleSeqs), as.character(sampleSeqs))
			if(exists("verbose")) cat(paste(length(unique(thisAlign$seq)),"/", thisAlign$nb,"unique seqs in original alignment, "))

			newAlign <- errorCorrect.alignment(thisAlign, correctThreshold)
			# need to repack this alignment

			if(exists("verbose")) cat(paste(length(unique(newAlign$seq)),"/", newAlign$nb,"unique seqs in new alignment, "))
			## DONE: repack the corrected alignment re-attribute the allele names. Update the markerSampleTable
			newCountTable <-  as.data.frame(table(unlist(newAlign$seq)),stringsAsFactors=FALSE)

			rawSeqCountTable <- data.frame(alignedVar=newCountTable$Var1, count=as.numeric(newCountTable$Freq),stringsAsFactors=FALSE)
			rawSeqCountTable$var <- gsub("-","",rawSeqCountTable$alignedVar)
			rawSeqCountTable <- rawSeqCountTable[order(rawSeqCountTable$count,decreasing=T),]
			#cat(paste("rawSeqCountTable rows:", nrow(rawSeqCountTable), "\n"))
			rawAlignFile.corrected <- paste("test", thisMarker, thisSample, "raw.align.corrected.fasta",sep=".")  #"localAlign.fasta"
			#write.fasta(as.list(newAlign$seq), newAlign$nam, file.out=rawAlignFile.corrected )
			write.fasta(as.list(rawSeqCountTable$alignedVar), rawSeqCountTable$var, file.out=rawAlignFile.corrected )
			rawAlignFile <- rawAlignFile.corrected
	}

	#v0.14# Align marker sequence using profile alignment (not done as part of base alignment so that errorCorrect can be run on rawAlignment).
	markerSeqFile <- paste(thisMarker, "markerSeq.fasta", sep=".")
	write.fasta(markerSeq, thisMarker, file.out=markerSeqFile )
	rawPlusMarkerFile <- paste("test", thisMarker, thisSample, "raw.marker.align.fasta",sep=".")  #"localAlign.fasta"
	muscleCommand <- paste(musclePath, "-profile -quiet -in1", rawAlignFile , "-in2", markerSeqFile  ,  "-out", rawPlusMarkerFile  )
	system(muscleCommand)

	#v0.14# Extract portion corresponding to reference. 
	rawAlignFile <- rawPlusMarkerFile 
	rawAlignment <- read.fasta(rawAlignFile, as.string=T)		# do not use read.alignment() - broken
	alignedMarkerSeq <- s2c(rawAlignment[[thisMarker]])
	subStart <- min(grep("-",alignedMarkerSeq ,invert=T))
	subEnd <- max(grep("-",alignedMarkerSeq ,invert=T))
	alignedSubSeqs <- lapply(rawAlignment, FUN=function(x)	substr(x[1], subStart, subEnd))
	subAlignFile <- paste("test", thisMarker, thisSample, "sub.align.fasta",sep=".")  #"localAlign.fasta"
	#subAlignFile <- paste(runPath, subAlignFileName , sep="/")
	write.fasta(alignedSubSeqs , names(alignedSubSeqs ), file.out=subAlignFile )

	alignedSubTable <- data.frame(var =  names(alignedSubSeqs ) , subSeq= as.character(unlist(alignedSubSeqs )))

	# Re-apply count of each seq.  There may be some duplicated subSeqs.
	combTable <- merge(rawSeqCountTable ,alignedSubTable , by="var", all.x=T)
	varCount <- by(combTable, as.character(combTable$subSeq), FUN=function(x) sum(x$count))
	varCountTable <- data.frame(alignedVar=names(varCount), count=as.numeric(varCount),stringsAsFactors=FALSE)	## added stringsAsFactors=FALSE to enable passing of aligned list
	varCountTable$var <- gsub("-","",varCountTable$alignedVar)
	varCountTable <- varCountTable[order(varCountTable$count,decreasing=T),]




	thisVarCount <- new("varCount", rawTotal=totalRaw,
			rawUniqueCount = rawVariantCount ,
			usedRawTotal = usedTotalRaw,
			usedRawUniqueCount = usedVarCount,
			subAlignTotal = sum(varCountTable$count),
			subAlignUniqueCount = nrow(varCountTable), 
			varCountTable = varCountTable)

	# Make unique list, summing counts where same seq found. (easier in table than list).  
	# ?TODO?
	#return(varCountTable)
	return(thisVarCount)

}







#mlgt <- function(object) attributes(object)
#setGeneric("mlgt")

#' Get variants for all markers/samples
#' 
#' \code{mlgt} Works through all pairs of markers and samples. Aligns variants and trims aligned variants to the marker sequence. Potential 'alleles' are assigned from the most common variants within each sample.
#'
#' Depends upon \code{\link{prepareMlgtRun}} having been run in the current directory to generate \option{designObject} of class \code{\link{mlgtDesign}}. 
#' The basic process for each marker/sample pair is to align all unique variants using MUSCLE and then extract the alignment portion aligned to the reference marker sequence, ignoring the rest.
#' The marker alignment is critical and \code{\link{mlgt}} has several options to optimise this alignment.
#' If the total number of reads is less than minTotalCount, then all variants are aligned. Otherwise, only the most abundant 30 unique variants are aligned.
#' Optionally, alignments are `error-correted' as per the separate function \code{\link{errorCorrect}}. Reads shorter than 'minLength' are filtered out.
#' 
#' @param designObject an object of class \code{\link{mlgtDesign}}
#' @param minTotalCount How many assigned sequences to allow before limiting the number of raw variants to allign.
#' @param maxVarsToAlign If total assigned sequences exceeds 'minTotalCount', then only the 'maxVarsToAlign' most abundant variants are used.
#' @param errorCorrect Use error correection on alignment of raw variants
#' @param correctThreshold Maximum proportion of raw reads at which (minor allele) bases and gaps are corrected.
#' @param minLength Reads below this length are excluded (they are very likely to be primer-dimers).
#'
#' @return an object of class \code{\link{mlgtResult}} containing all variants and their counts, a summary table (all markers) and one summary table per marker.
#' @seealso \code{\link{prepareMlgtRun}}
#'
#' @export
#' @docType methods
#' @rdname mlgt-methods
#' @aliases mlgt.mlgtDesign
mlgt <- function(designObject, maxVarsToAlign=30, minTotalCount=500, errorCorrect=FALSE,correctThreshold=0.01, minLength=70) attributes(designObject)
setGeneric("mlgt")


mlgt.mlgtDesign <- function(designObject, maxVarsToAlign=30, minTotalCount=500, errorCorrect=FALSE,correctThreshold=0.01, minLength=70)  {
	topHits <- getTopBlastHits("blastOut.markers.tab")
	topHits$strand <- ifelse(topHits$s.end > topHits$s.start, 1,2)
	fMarkerMap <- split(as.character(topHits$query[topHits$strand==1]), topHits$subject[topHits$strand==1])
	rMarkerMap <- split(as.character(topHits$query[topHits$strand==2]), topHits$subject[topHits$strand==2])

	## NEED TO MAKE SAMPLEMAP WITH HITS TO MID IN BOTH FORWARD AND REVERSE STRANDS like marker hits are split.
	## Requires retention of 2 blast hits per sequence.
	topSampleHits <- read.delim("blastOut.rTags.tab", header=F)
	names(topSampleHits ) <- c("query", "subject", "percentId", "aliLength", "mismatches", "gapOpenings", "q.start","q.end", "s.start","s.end", "p_value", "e_value")
	topSampleHits$strand <- ifelse(topSampleHits$s.end > topSampleHits$s.start, 1,2)
	fSampleMap <- split(as.character(topSampleHits$query[topSampleHits$strand==1]), topSampleHits$subject[topSampleHits$strand==1])
	rSampleMap <- split(as.character(topSampleHits$query[topSampleHits$strand==2]), topSampleHits$subject[topSampleHits$strand==2])
	# combind sampleMaps to give sequences with MIDs in both orientations.
	pairedSampleMap <- lapply(names(fSampleMap), FUN=function(x) intersect(fSampleMap[[x]], rSampleMap[[x]]))
	names(pairedSampleMap) <- names(fSampleMap)

	##########ITERATIONS

	markerSampleList <- list()
	runSummaryTable <- data.frame()
	alleleDb <- list()
	varCountTableList <- list()

	if(errorCorrect)  cat(paste("Using error correction at", correctThreshold, "\n"))

	for(thisMarker in names(designObject@markers)) {
	#for(thisMarker in names(markerMap)) {
	#for(thisMarker in names(markerMap)[1:2]) {	# temp to finish off

	cat(paste(thisMarker,"\n"))
	#cat("New Version\n")
	#thisMarker <- "DQA1_E2"

	## might need to combine all these to return a single item.
	summaryList <- list()
	summaryTable <- data.frame()
	markerSequenceCount <- list("noSeq"=0)		#  BUG? requires some data otherwise won't sum properly with localSequenceCount.
	alleleList <- list() 
	variantList <- list()
	alleleCount <- 1
	markerSeq <- unlist(getSequence(designObject@markers[[thisMarker]],as.string=T))
	varCountTableList[[thisMarker]] <- data.frame()
	
	for(thisSample in designObject@samples) {
		#for(thisSample in names(pairedSampleMap)[1:4]) {
			#print(thisSample)

			testPairSeqList <- intersect(pairedSampleMap[[thisSample]],union(fMarkerMap[[thisMarker]], rMarkerMap[[thisMarker]]))
			#testPairSeqList <- intersect(pairedSampleMap[[thisSample]], markerMap[[thisMarker]])
			#testPairSeqList <- intersect(sampleMap[[thisSample]], markerMap[[thisMarker]])


		seqTable <- data.frame()
		localAlleleNames <- c("NA","NA","NA")
		localAlleleFreqs <- c(0,0,0)

		## go through all seq's mapped to this marker/sample pair.
		## extract the corresponding sequence delimited by the top blast hits on the primers.  IS THIS THE BEST WAY?
		##		Simple improvement: minimum blast hit length to primer to keep. 

		## internal Function

		recordNoSeqs <- function(summaryTable)  {		# to record no seqs before skipping out. 
				summaryRow <- data.frame(marker=thisMarker, sample=thisSample, 
					rawTotal=0, rawVars=0,
					usedTotal=0, usedVars=0,
					numbSeqs=0,numbVars=0,
					varName.1="NA", varFreq.1= 0,
					varName.2="NA", varFreq.2= 0,
					varName.3="NA", varFreq.3= 0)
				summaryTable <- rbind(summaryTable, summaryRow)
				return(summaryTable)
		}



		if(length(testPairSeqList) < 1) {
			#summaryList[[thisMarker]][[thisSample]] <- NA	
			summaryTable  <- recordNoSeqs(summaryTable)
			next ;	# skip to next sample
		} 

		## Ver. 0.14 - edited getSubSeqsTable to return object of class 'varCount' , which includes the required table.
		# seqTable <- getSubSeqsTable(thisMarker, thisSample, pairedSampleMap, fMarkerMap,rMarkerMap, markerSeq)
		thisVarCount <-  getSubSeqsTable(thisMarker, thisSample, pairedSampleMap, fMarkerMap,rMarkerMap, markerSeq, 
					maxVarsToAlign=maxVarsToAlign, minTotalCount=minTotalCount, errorCorrect=errorCorrect, 
					correctThreshold=correctThreshold, minLength=minLength)
		seqTable <- thisVarCount@varCountTable		

		#cat(paste("Raw total:",thisVarCount@rawTotal,"\n"))

		# if no sequences returned, nothing to process. 
		if(nrow(seqTable) < 1 )  {
			summaryTable  <- recordNoSeqs(summaryTable)
			#summaryList[[thisMarker]][[thisSample]] <- NA	
			next ;		# go to next sample.
		}

		# store sequence and count of sequence as alignedVar (needed for alignment report)
		varCountTableList[[thisMarker]][seqTable$alignedVar,thisSample] <- seqTable$count
		#varCountTableList[[thisMarker]][seqTable$var,thisSample] <- seqTable$count

		#localSequenceMap <- split(seqTable[,2], seqTable[,1])

		#localSequenceCount <- lapply(localSequenceMap , length)  # list named by sequence with counts.
		#localSequenceCount <- localSequenceCount[order(as.numeric(localSequenceCount), decreasing=T)]

		## test if variants are novel. 
		## Give allele names?  
		## Do with first three for now. 


		alToRecord <- min(3,nrow(seqTable))
		if(alToRecord > 0)  {
			for (a in 1:alToRecord )  {
				if(is.null(variantList[[seqTable$var[a]]]))  {    	# novel
					alleleName <- paste(thisMarker, alleleCount,sep=".")	
					variantList[[seqTable$var[a]]] <- alleleName
					localAlleleNames[a] <- alleleName 
					localAlleleFreqs[a] <- seqTable$count[a]
					alleleCount <- alleleCount + 1
				} else  {										# pre-existing alllele
					localAlleleNames[a] <- variantList[[seqTable$var[a]]]
					localAlleleFreqs[a] <- seqTable$count[a]		
				}
			}
		}


		# sequence correction?  


		# compile stats

		if(nrow(seqTable) >0 )  {	# cannot allow assignment from empty list as messes up class of list for remaining iterations
			summaryList[[thisMarker]] <- list()
			summaryList[[thisMarker]][[thisSample]] <- seqTable
		}

		summaryRow <- data.frame(marker=thisMarker, sample=thisSample, 
					rawTotal=thisVarCount@rawTotal, rawVars=thisVarCount@rawUniqueCount,
					usedTotal=thisVarCount@usedRawTotal, usedVars=thisVarCount@usedRawUniqueCount,
					numbSeqs=sum(seqTable$count),numbVars=nrow(seqTable),
					varName.1=localAlleleNames[1], varFreq.1= localAlleleFreqs[1],
					varName.2=localAlleleNames[2], varFreq.2= localAlleleFreqs[2],
					varName.3=localAlleleNames[3], varFreq.3= localAlleleFreqs[3])
		summaryTable <- rbind(summaryTable, summaryRow)

		#sequence count across samples? 
		# need to sum from summaryTable or from summaryList.
		#markerSequenceCount <- 
		#as.list(colSums(merge(m, n, all = TRUE), na.rm = TRUE))  # not working
		localSequenceCount <- as.list(seqTable$count)
		names(localSequenceCount) <- seqTable$var
		markerSequenceCount   <- as.list(colSums(merge(markerSequenceCount  , localSequenceCount,  all = TRUE), na.rm = TRUE))
		# might need to instantiate the markerSequenceCount if empty. 



	}  # end of sample loop

	markerSampleList[[thisMarker]] <- summaryTable
	## DONE: min(nchar(names(variantList))) throws warning when no variants in list. (Inf/-Inf)
	minVarLength <- ifelse(length(markerSequenceCount) < 1, NA, min(nchar(names(markerSequenceCount))) )
	maxVarLength <- ifelse(length(markerSequenceCount) < 1, NA, max(nchar(names(markerSequenceCount))) )
	minAleLength <- ifelse(length(variantList) < 1, NA, min(nchar(names(variantList))) )
	maxAleLength <- ifelse(length(variantList) < 1, NA, max(nchar(names(variantList))) )

	runSummaryRow <- data.frame(marker=thisMarker, assignedSeqs=sum(summaryTable$numbSeqs), assignedVariants=sum(summaryTable$numbVars), 
					minVariantLength=minVarLength, 
					maxVariantLength=maxVarLength,
					minAlleleLength=minAleLength, maxAlleleLength=maxAleLength )
	runSummaryTable <- rbind(runSummaryTable, runSummaryRow)
	if(length(variantList) > 0)  {
		# This line replaced. Not entirely tested the repurcussions. e.g. makeVarAlleleMap()?
		#alleleDb[[thisMarker]] <- variantList   # LATER: separate lists for alleles and variants? 
		#alleleDb[[thisMarker]] <- list(reference=as.SeqFastadna(markerSeq, thisMarker), alleleMap=variantList, inputAlleleCount = length(unlist(variantList)), uniqueSubAlleleCount=length(variantList))
		alleleDb[[thisMarker]] <- new("variantMap", reference=as.SeqFastadna(markerSeq, thisMarker), 
							variantSource=paste(designObject@projectName, designObject@runName,sep="."),
							variantMap=variantList, inputVariantCount = length(unlist(variantList)), uniqueSubVariantCount=length(variantList))
	}

}  # end of marker loop
	
	localMlgtResult <- new("mlgtResult", designObject,  runSummaryTable=runSummaryTable , alleleDb=alleleDb, markerSampleList=markerSampleList,
						varCountTables=varCountTableList)
	return(localMlgtResult)

}  # end of mlgt function

#' @rdname mlgt-methods
#' @aliases mlgt,mlgtDesign-method
setMethod("mlgt",signature(designObject="mlgtDesign", maxVarsToAlign="ANY", minTotalCount="ANY", errorCorrect="ANY",correctThreshold="ANY", minLength="ANY"), definition=mlgt.mlgtDesign)

#INTERNAL. Does this need documenting?
# Create a local BLAST db 
# 
# This internal utility uses \code{system(fomatdb)} to create a local BLAST database
#
# Requires the NCBI program formatdb to be installed and the environment variable formatdbPath to be set.
#
# @param inputFastaFile
# @param formatdbPath
# @param blastdbName 
# @param indexDb Whether to generate an index for the db. Useful for large db of sequences, which are later to be extracted with fastacommand
#
#
setUpBlastDb <- function(inputFastaFile, formatdbPath, blastdbName, indexDb="F")  {
	formatdbCommand <- paste(formatdbPath, "-i",  inputFastaFile ,  "-p F -o", indexDb ,"-n", blastdbName)
	#cat(paste(formatdbCommand,"\n"))
	system(formatdbCommand)
}


quoteIfSpaces <- function(pathString)  {
	if(length(grep(" ",pathString, fixed=T))  > 0 ) {
		  pathString <- shQuote(pathString)	
	}
	return(pathString)
}

copyIfSpaces <- function(pathString)  {
	if(length(grep(" ",pathString, fixed=T))  > 0 ) {
		newName <- make.names(basename(pathString))
		file.copy(pathString,newName , overwrite=TRUE)
		cat(paste("BLAST and MUSCLE cannot cope with whitespace in filenames.\nCopying",pathString,"to",newName,"\n"))
		return(newName)
	} else {
		return(pathString)	
	}
}


#prepareMlgtRun <- function(object) attributes(object)
#prepareMlgtRun <- function(designObject) attributes(designObject)
#' Prepare to run mlgt
#' 
#' Required before \code{\link{mlgt}} is used. Create BLAST databases and assign sequences using BLAST. 
#'
#' This important function stores all the information about the analysis run AND populates the working directory 
#' with multiple local Blast databases, which are later required by \code{\link{mlgt}}. 
#' Once \code{prepareMlgtRun} has been run,  \code{\link{mlgt}} can be run aswell as 
#' \code{\link{printBlastResultGraphs}} and \code{\link{inspectBlastResults}}.
#' 
#' @param designObject Only used internally.
#' @param projectName In which project does this run belong
#' @param runName Which run was this. An identifier for the sequnece run
#' @param markers	A \emph{list} of named sequences. 
#' @param samples A vector of sample names 
#' @param fTags A vector of named sequence of MIDs used to barcode samples at the 5' end.
#' @param rTags A vector of named sequence of MIDs used to barcode samples at the 3' end.
#' @param inputFastaFile The name of the file containing sequences. Currently only fasta format is supported. It is up to you to pre-filter the sequences.
#' @param overwrite Should files in the current directory be overwritten? c("prompt", "yes", "no")
#'
#' @return An object of class \code{\link{mlgtDesign}} is returned. Also, several BLAST dbs and sets of BLAST results are created in the working directory.
#'	These are essential for \code{\link{mlgt}} to run.
#'
#' @rdname prepareMlgtRun-methods
#' @export
#' @aliases prepareMlgtRun.listDesign,prepareMlgtRun.mlgtDesign 
#' @seealso \code{\link{printBlastResultGraphs}} and \code{\link{inspectBlastResults}} can only be run AFTER \code{prepareMlgtRun}.
prepareMlgtRun <- function(designObject,projectName,runName, samples, markers,fTags,rTags, inputFastaFile, overwrite) attributes(designObject)

setGeneric("prepareMlgtRun")

prepareMlgtRun.listDesign <- function(projectName,runName, samples, markers,fTags,rTags, inputFastaFile,overwrite="prompt")  {
	# test inputFastaFile for spaces 
	
	inputFastaFile <- copyIfSpaces(inputFastaFile)
	#if(length(grep(" ",inputFastaFile, fixed=T))  > 0 ) stop(paste("mlgt cannot handle path names with whitespace. Sorry.\n ->>",inputFastaFile,"\n"))
	##inputFastaFile <- quoteIfSpaces(inputFastaFile)

	designObject <- new("mlgtDesign", projectName=projectName, runName=runName, 
				samples=samples, markers=markers ,
				fTags=fTags, rTags=rTags, inputFastaFile=inputFastaFile)

	designObject <- prepareMlgtRun(designObject,overwrite=overwrite)
	
}

#' @rdname prepareMlgtRun-methods
#' @aliases prepareMlgtRun,missing,character,character,character,list,list,list,character,character-method
setMethod("prepareMlgtRun",
	signature(designObject="missing",projectName="character", runName="character", samples="character",markers="list", 
	fTags="list", rTags="list", inputFastaFile="character", overwrite="character"), 
	definition=prepareMlgtRun.listDesign)




prepareMlgtRun.mlgtDesign <- function(designObject, overwrite="prompt")  {
	cat(paste(designObject@projectName,"\n"))

	#check pathNames of auxillary programs.
	pathNames <- c("BLASTALL_PATH","FORMATDB_PATH","FASTACMD_PATH","MUSCLE_PATH")
	for(thisPath in pathNames)  {
		if(nchar(Sys.getenv(thisPath)) < 1) {
			stop(paste(thisPath,"has not been set!"))
		}
		# shellQuote any paths containing spaces.
		# Can't use variables in Sys.setenv
		#if(length(grep(" ",Sys.getenv(thisPath), fixed=T))  > 0 )  Sys.setenv(thisPath= shQuote(thisPath))
	}

	formatdbPath <- Sys.getenv("FORMATDB_PATH")
	fastacmdPath <- Sys.getenv("FASTACMD_PATH")
	blastAllPath <- Sys.getenv("BLASTALL_PATH")

	# shellQuote any paths containing spaces.
	if(length(grep(" ",formatdbPath, fixed=T))  > 0 ) formatdbPath <- shQuote(formatdbPath)
	if(length(grep(" ",fastacmdPath, fixed=T))  > 0 ) fastacmdPath <- shQuote(fastacmdPath)
	if(length(grep(" ",blastAllPath, fixed=T))  > 0 ) blastAllPath <- shQuote(blastAllPath)



	# runPath <- getwd()	# could change this later
	# if(length(grep(" ",runPath , fixed=T))  > 0 )  {
	# 	stop("mlgt cannot yet cope with path names containing spaces. Please move to a folder with no spaces in path. Sorry. \n")
	# }
	# fTagsFastaFile <- paste(runPath,"fTags.fasta", sep="/")
	# rTagsFastaFile <- paste(runPath,"rTags.fasta", sep="/")
	# markersFastaFile <- paste(runPath,"markers.fasta", sep="/")
	# rTagsBlastOutFile <- paste(runPath,"blastOut.rTags.tab", sep="/")
	# fTagsBlastOutFile <- paste(runPath,"blastOut.fTags.tab", sep="/")
	# blastOutFileName <- "blastOut.markers.tab"
	# markerBlastOutFile <- paste(runPath, blastOutFileName, sep="/")

	# removed concatenation of runPath because of issue with space character. All output must now go to working directory
	fTagsFastaFile <- "fTags.fasta"
	rTagsFastaFile <- "rTags.fasta"
	markersFastaFile <- "markers.fasta"
	rTagsBlastOutFile <- "blastOut.rTags.tab"
	fTagsBlastOutFile <- "blastOut.fTags.tab"
	blastOutFileName <- "blastOut.markers.tab"
	markerBlastOutFile <- blastOutFileName	#redundant?
	existingFiles <- (file.exists(fTagsFastaFile ) | file.exists(rTagsFastaFile ) | file.exists(markersFastaFile ) | 
					file.exists(rTagsBlastOutFile) | file.exists(fTagsBlastOutFile) | file.exists(markerBlastOutFile))
	if(existingFiles)  {
		overwrite <- tolower(overwrite)
		# set up new directories if required.
		if(overwrite == "prompt") overwrite <- readline("This folder already contains mlgt run files, do you want to write over this data? (yes/no)")

		if(!(overwrite=="y" | overwrite=="n" | overwrite=="yes" | overwrite=="no")) {stop(paste("Unrecognised value for overwrite:", overwrite))}
		overWriteBaseData <- switch(overwrite,
						"yes"=TRUE,
						"y"=TRUE,
						"no"=FALSE,
						"n"=FALSE)
		if(!overWriteBaseData) {stop("This folder already contains mlgt run files. Exiting")}
	}					

	#runPath <- getwd()
	# set up blast DBs
	cat("Setting up BLAST DBs...\n")

	write.fasta(designObject@fTags, names(designObject@fTags), file.out=fTagsFastaFile)
	setUpBlastDb(fTagsFastaFile , formatdbPath, blastdbName="fTags")		# default is no index

	write.fasta(designObject@rTags, names(designObject@rTags), file.out=rTagsFastaFile )
	setUpBlastDb(rTagsFastaFile , formatdbPath, blastdbName="rTags")		# default is no index

	write.fasta(designObject@markers, names(designObject@markers), file.out=markersFastaFile )
	setUpBlastDb(markersFastaFile , formatdbPath, blastdbName="markerSeqs", indexDb="T") 

	## input as blast DB with index (useful for fast sub-sequence retrieval?)
	inputFastaFile <- designObject@inputFastaFile		#
	setUpBlastDb(inputFastaFile , formatdbPath, blastdbName="inputSeqs", indexDb="T") 

	# run preliminary blast
	cat("Running BLAST searches...\n")
	
	rMinTagSize <- 10	# limit blast hits to perfect matches.		### TODO: SET GLOBALLY
	dbName <- "rTags"

	blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", inputFastaFile , "-W", rMinTagSize  ,"-m 8 -b 2 -S 3 -o", rTagsBlastOutFile )
	system(blastCommand )

	fMinTagSize <- 10	# limit blast hits to perfect matches.
	dbName <- "fTags"

	blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", inputFastaFile , "-W", fMinTagSize ,"-m 8 -b 2 -S 3 -o", fTagsBlastOutFile )
	system(blastCommand )

	## blast against markers
	dbName <- "markerSeqs"

	blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", inputFastaFile , "-W", 11 , "-m 8 -b 1 -o", markerBlastOutFile )		
	system(blastCommand )

	designObject@markerBlastResults <- markerBlastOutFile
	designObject@fTagBlastResults <- fTagsBlastOutFile 
	designObject@rTagBlastResults <- rTagsBlastOutFile 

	return(designObject)

}	

##setMethod("prepareMlgtRun","mlgtDesign", definition=prepareMlgtRun.mlgtDesign)
#setMethod("prepareMlgtRun",signature(designObject="mlgtDesign"), definition=prepareMlgtRun.mlgtDesign)


#TODO: test if signature with overwrite="ANY" will work. First attempt, no.
#' @rdname prepareMlgtRun-methods
#' @aliases prepareMlgtRun,mlgtDesign,missing,missing,missing,missing,missing,missing,missing,character-method
setMethod("prepareMlgtRun",
	signature(designObject="mlgtDesign",projectName="missing", runName="missing", samples="missing",markers="missing", 
	fTags="missing", rTags="missing", inputFastaFile="missing", overwrite="character"),
	 definition=prepareMlgtRun.mlgtDesign)


#setGeneric("prepareMlgtRun","mlgtDesign", definition=prepareMlgtRun.mlgtDesign)

######################### genotyping & allele matching

#' Results from \code{\link{callGenotypes}}
#'
#' An S4 class containing a table and parameter values returned by \code{\link{callGenotypes}}
#'
#' \describe{
#'   \item{projectName}{In which project does this run belong}
#'   \item{runName}{Which run was this. An identifier for the sequence run}
#'   \item{marker}{Which marker was this.} 
#'   \item{genotypeTable}{A data frame with variant counts, statistics, genotype calls and, optionally, allele names.}
#'   \item{callMethod}{Which method was used to call genotypes}
#'   \item{callParameters}{a named list containing parameter values used in the call}
#'   \item{mappedToAlleles}{TRUE/FALSE whether an attempt was made to map the variants to a db on known alleles.}
#'   \item{alleleDbName}{A list of objects of class \code{\link{variantMap}}. Contains all variants returned by \code{\link{mlgt}}}
#'  }
#'
#' @seealso \code{\link{callGenotypes}}, \code{\link{writeGenotypeCallsToFile}}
setClass("genotypeCall", 
	representation(
		projectName="character", 
		runName="character",
		marker="character",
		genotypeTable="data.frame",
		callMethod="character",
		callParameters="list",
		mappedToAlleles="logical",
		alleleDbName="character"
	)
)

#INTERNAL. Does this need documentation?
makeVarAlleleMap <- function(allele.variantMap, variant.variantMap)  {
			varAlleleMap <- data.frame()
			# to use stopifnot, need to ensure no empty maps are passed. Currently this is happening so taking out this check.
			#stopifnot(allele.variantMap@reference == variant.variantMap@reference)
			knownAlleleTable <- data.frame(alleleSeq=names(allele.variantMap@variantMap), knownAlleles=as.character(allele.variantMap@variantMap))


			dataAlleleTable <-  data.frame(alleleSeq=names(variant.variantMap@variantMap), varNames=as.character(variant.variantMap@variantMap))
			
			varAlleleMap <- merge(knownAlleleTable, dataAlleleTable , by="alleleSeq")
}

# deprecated and/or defunct
makeVarAlleleMap.list <- function(alleleDb, varDb,alleleMarkers=names(alleleDb),varMarkers=names(varDb))  {
			knownAlleleTable <- data.frame()
			for(thisMarker in alleleMarkers)  {
				knownAlleleTable <- rbind(knownAlleleTable , data.frame(alleleSeq=names(alleleDb[[thisMarker]]@alleleMap), knownAlleles=as.character(alleleDb[[thisMarker]]@alleleMap)))
			}

			dataAlleleTable <- data.frame()
			for(thisMarker in varMarkers)  {
				# this first line is what it SHOULD be like, once mlgtResult is updated to match the alleleDb format. Need new class: alleleDb
				dataAlleleTable <- rbind(dataAlleleTable , data.frame(alleleSeq=names(varDb[[thisMarker]]@alleleMap), varNames=as.character(varDb[[thisMarker]]@alleleMap)))
				#dataAlleleTable <- rbind(dataAlleleTable , data.frame(alleleSeq=names(varDb[[thisMarker]]), varNames=as.character(varDb[[thisMarker]])))
			}
			
			varAlleleMap <- merge(knownAlleleTable, dataAlleleTable , by="alleleSeq")
}


vectorOrRepeat <- function(paramValue, requiredLength) {
	# Used by callGenotypes.mgltResult() to equalise lengths of parameter vectors when one has length > 1
	if(length(paramValue) == requiredLength) {
		return(paramValue)
	} else {
		if(length(paramValue) == 1) {
			return(rep(paramValue, requiredLength))
		} else {
			stop(paste("Parameter has length",length(paramValue) ,"but does not match required length", requiredLength))	
		}
	}

}

# replacement for vectorOrRepeat for cases where custom methods are allowed in callGenotypes().
# Returns a list of lists each identical except where user specified a vector of parameter values.
makeBigParamList <- function(..., markerCount)  {
	params <- list(...)
	paramBigList <- list()
	for(i in 1:markerCount) {
		paramBigList[[i]] <- params 
	}
	lengthList <- lapply(params , length)
	vecParams <- which(lengthList > 1)
	if(length(vecParams) > 0) {
		vecNames <- names(params)[vecParams]
		#cat(vecParams)
		if(any(lengthList[vecParams] != markerCount)) { stop("Vector supplied that does not match length of markerList") }
		for(i in 1:length(vecParams)) {
			for(j in 1:markerCount) {
				paramBigList[[j]][[vecNames[i]]] <- params[[vecNames[i]]][j] 
			}
		}
	}
	return(paramBigList)
}


# Used for approximate (BLAST) variant to allele matching.
# returns a standard blast table with single best hit per query.
makeVarAlleleBlastMap <- function(allele.variantMap, variant.variantMap)  {
	pathNames <- c("BLASTALL_PATH","FORMATDB_PATH","FASTACMD_PATH","MUSCLE_PATH")
	for(thisPath in pathNames)  {
		if(nchar(Sys.getenv(thisPath)) < 1) {
			stop(paste(thisPath,"has not been set!"))
		}
		# shellQuote any paths containing spaces.
		# Can't use variables in Sys.setenv
		#if(length(grep(" ",Sys.getenv(thisPath), fixed=T))  > 0 )  Sys.setenv(thisPath= shQuote(thisPath))
	}

	formatdbPath <- Sys.getenv("FORMATDB_PATH")
	fastacmdPath <- Sys.getenv("FASTACMD_PATH")
	blastAllPath <- Sys.getenv("BLASTALL_PATH")

	# shellQuote any paths containing spaces.
	if(length(grep(" ",formatdbPath, fixed=T))  > 0 ) formatdbPath <- shQuote(formatdbPath)
	if(length(grep(" ",fastacmdPath, fixed=T))  > 0 ) fastacmdPath <- shQuote(fastacmdPath)
	if(length(grep(" ",blastAllPath, fixed=T))  > 0 ) blastAllPath <- shQuote(blastAllPath)

	thisMarker <-  getName(allele.variantMap@reference)
	cat(thisMarker)
	#export known variants to fasta (may need to sort out names)
	knownAsFastaFile <- paste(thisMarker, "knownAlleles.fasta", sep=".")	
	# some sub-alleles may have zero length (they didn't align at all with marker)
	thisVariantMap <- allele.variantMap@variantMap
	thisVariantMap <- thisVariantMap[nchar(names(thisVariantMap)) > 0]
	#write the known alleles to a fasta file. The concatentated names must be truncated to XX chars for blast to run.
	write.fasta(as.list(names(thisVariantMap)),substring(as.character(thisVariantMap),1,60), file.out=knownAsFastaFile )
	#create blast DB of known variants.
	dbName <- paste(thisMarker,"knownAlleles",sep=".")
	setUpBlastDb(knownAsFastaFile , formatdbPath, blastdbName=dbName )
	#export new variants to fasta
	newVarsAsFastaFile <- paste(thisMarker, "newVars.fasta", sep=".")
	write.fasta(as.list(names(variant.variantMap@variantMap)),as.character(variant.variantMap@variantMap), file.out=newVarsAsFastaFile )
	#blast new variants against known DB
	newVarBlastResultFile <- paste(thisMarker, "blastOut.newVars.tab", sep=".")
	blastCommand <- paste(blastAllPath , "-p blastn -d", dbName , "-i", newVarsAsFastaFile , "-W", 11 , "-m 8 -b 1 -o", newVarBlastResultFile)
	system(blastCommand)
	#open and retrieve the results.
	topHits <- getTopBlastHits(newVarBlastResultFile)
	#return as a useful lookuptable
}


## call genotypes on a table of variant counts. Can select specific markers/samples to return. 
#' Default internal methods for \code{\link{callGenotypes}} 
#'
#' This is the default method to call genotypes from a table of variant counts.
#' Methods:-
#' \describe{
#' 	\item{`callGenotypes.default'}{Three sequential steps for each marker/sample pair: 
#' 		\enumerate{
#'			\item {if the number of reads is less than \code{minTotalReads} the genotype is \emph{`tooFewReads'} }
#'			\item {if the difference between the sum of counts of the top two variants and the count of the third most variant, expressed as proportion of total, is less than \code{minDiffToVarThree}, OR the third most abundant variant accounts for more than maxPropVarThree (default=0.1) of the reads, then the genotype is \emph{`complexVars'}}
#'			\item {if the difference between the counts of top two variants, expressed as a proportion of the total, is greater than or equal to \code{minPropDiffHomHetThreshold}, then the genotype is \emph{HOMOZYGOTE}. Otherwise it is \emph{HETEROZYGOTE}. }  
#'		}
#' 	}
#' }
#'
#' @param table The table of sequence counts as in the markerSampleTable of an mlgtResult object.
#' @param minTotalReads Minimum number of reads before attempting to call genotypes
#' @param minDiffToVarThree Difference between sum of counts of top two variants and the count of the third most frequent variant, expressed as proportion of total. 
#' @param minPropDiffHomHetThreshold Difference between counts of top two variants. One way to distinguish HOMOZYGOTES and HETEROZYGOTES.
#' @param maxPropVarThree Also call as 'complexVars' if the third variant accounts for more than this proportion of used reads (default=0.1)
#' @return A data.frame identical to those in markerSampleList but with additional columns giving parameter values, 
#' and a 'status' column giving the genotype status.
callGenotypes.default <- function(table,  minTotalReads=50, minDiffToVarThree=0.4,
					minPropDiffHomHetThreshold=0.3, maxPropVarThree=0.1) {
	
	#table$genotype
	table$status <- "notCalled"
	enoughReads <- table$numbSeqs >= minTotalReads
	table$status[!enoughReads] <- "tooFewReads"
	
	# difference between sum of vars 1 + 2 and var 3, as proportion of total
	# > 0.5 is good. <0.3 not good. 0.4 as cut-off for now? 
	table$diffToVarThree <- with(table, ((varFreq.1+varFreq.2)-varFreq.3)/numbSeqs)
	table$propThree <- with(table, (varFreq.3)/numbSeqs)

	#distinctVars <- 	with(table, diffToVarThree  >= minDiffToVarThree)
	distinctVars <- with(table, (diffToVarThree  >= minDiffToVarThree)  & propThree < maxPropVarThree)
	table$status[enoughReads & !distinctVars] <- "complexVars"

	# difference between var 1 and var2 as proportion of total
	# homozygote: >0.3, often > 0.5
	# heterozygote: <0.25, often < 0.1
	table$propDiffHomHet <- with(table, ((varFreq.1-varFreq.2)/numbSeqs))
	eligible <- (enoughReads & distinctVars) 
	table$status[eligible] <- ifelse((table$propDiffHomHet[eligible]  >= minPropDiffHomHetThreshold), "HOMOZYGOTE","HETEROZYGOTE")

	#table$homozygote <- 	with(table, ((varFreq.1-varFreq.2)/numbSeqs) >= minPropDiffHomHetThreshold)

	return(table)
}

## example custom function
callGenotypes.custom <- function(table) {
	table$status <- "notCalled"
	return(table)
}





## To make this customisable, would need to include '...' as first argument so that users can add in variables. T
## This function would then pass those variables to the user supplied function 
## Some rules will have to be made on required input and output. 
## e.g. result should still be list of class "genotypeCall" but used could specify what the table can contain.
## and the user supplied parameters could be stored within the callParameters slot.
## IMPORTANT, would also need to sort out the generic set-up given below (probably no longer needed). 
# generic method for callGenotypes.
#' Make genotype calls
#' 
#' Apply a genotype call method to a table or list of tables of variant data such as the \code{markerSampleList} table of an \code{\link{mlgtResult}}.
#'
#' After \code{\link{mlgt}} has generated tables of the most common variants assigned in each marker/sample pair, an attempt can be made to call genotypes.
#' This is kept separate because users might want to try different calling methods and have the option to map to a known set of alleles. Currently, only 
#' one method is implemented (\emph{`custom'}). See \code{\link{callGenotypes.default}}.
#' This function also includes the option to map variants to a list of known alleles created using \code{\link{createKnownAlleleList}}. The basic method makes only perfect matches but a secondary method can be triggered (approxMatching=TRUE) to find the allele with the greatest similarity using a local BLAST search.
#' 
#' 
#' @param resultObject An object of class \code{\link{mlgtResult}}, as returned by \code{\link{mlgt}}
#' @param alleleDb A list of \code{\link{variantMap}} objects derived from known alleles. As made by 
#' \code{\link{createKnownAlleleList}}
#' @param method How to call genotypes. Currently only "callGenotypes.default" is implemented. Users can define
#' their own methods as R functions (see the vignette).
#' @param markerList For which of the markers do you want to call genotypes (default is all)?
#' @param sampleList For which of the samples do you want to call genotypes (default is all)?
#' @param mapAlleles FALSE/TRUE. Whether to map variants to db \option{alleleDb} of known alleles. 
#' @param approxMatching If TRUE, a BLAST search is also performed to find matches (slower). Additional columns are added to the genoytpeTable
#' @param ... Other parameter values will be passed to custom methods such as \code{\link{callGenotypes.default}}
#'
#' @return list of call results including the call parameters and a table of calls (class \code{\link{genotypeCall}}). If an mlgtResult object was supplied then a list of \code{\link{genotypeCall}} objects will be returned, each named by marker.
#'
#' @export
#' @docType methods
#' @aliases callGenotypes.mlgtResult
#' @examples \dontrun{ 
#' data("mlgtResult", package="mlgt")
#' my.mlgt.Result
#' # the default method
#' my.genoytpes <- callGenotypes(my.mlgt.Result)
#' # using a custom method
#' callGenotypes.custom <- function(table, maxPropUniqueVars=0.5) {
#' 	table$status <- "notCalled"
#' 	table$propUniqueVars <- table$numbVar/table$numbSeq
#' 	table$status <- ifelse(table$propUniqueVars <= maxPropUniqueVars,"good", "bad")
#' 	return(table)
#' }
#' my.custom.Genotypes <- callGenotypes(my.mlgt.Result, method="callGenotypes.custom")
#' }
callGenotypes <- function(resultObject,  method="callGenotypes.default",  
					markerList=names(resultObject@markers),
					sampleList=resultObject@samples, mapAlleles=FALSE, alleleDb=NULL, approxMatching=FALSE, ...	) {

	## FM requested marker specific parameters.
	## test for vectors in any of the calling parameters. 	
	## if all are single, apply genotyping to one large table of all results,	
	## if any are vectors, need to apply genotyping to each marker separately. 
	## Should mark if different thresholds used. 
	## need to test that the threshold vector matches the markerlist length.

	runTable <- data.frame()
	genotypeTable <- data.frame()
	callResults <- list()

#	if(length(minTotalReads) > 1 | length(minDiffToVarThree) > 1 | length(minPropDiffHomHetThreshold) > 1 )  {
		# find parameters as vectors and test the lengths. Set those which are only 1 to be length of markerList.
		#minTotalReads <- vectorOrRepeat(minTotalReads, length(markerList))
		#minDiffToVarThree <- vectorOrRepeat(minDiffToVarThree, length(markerList))
		#minPropDiffHomHetThreshold<- vectorOrRepeat(minPropDiffHomHetThreshold, length(markerList))		

		longParamList <- makeBigParamList(..., markerCount=length(markerList))
		
		# multiple parameter values set. Genotype each marker separately
		for(i in 1:length(markerList))  {	
			thisMarker <- markerList[i]	
			subTable <- resultObject@markerSampleList[[thisMarker]]
			subTable <- subTable[match(sampleList,subTable$sample),]	
			if(nrow(subTable) < 1)  {
				# do nothing with this marker
				warning(paste("No data for:",thisMarker))
			} else {
				thisParamList <- longParamList[[i]]
				thisParamList[['table']] <- subTable
				genotypeTable <- do.call(method, thisParamList )
				
				#genotypeTable <- callGenotypes.table(subTable , alleleDb=alleleDb, method=method,minTotalReads=minTotalReads[i], 
				#	minDiffToVarThree=minDiffToVarThree[i],
				#	minPropDiffHomHetThreshold=minPropDiffHomHetThreshold[i], mapAlleles=mapAlleles)
				#genotypeTable <- rbind(genotypeTable , subGenoTable)

				if(mapAlleles) {
					# 2-step strategy. 
					# 1. Try perfect matching first. Quick and perfect.
					# Keep track of which alleles were matched this way
					# 2. Match with BLAST (build blast DB of known alleles and blast all new variants (once) against this.
					# Record best hit and quality of best hit (percentID, percentlength)
					if(is.null(alleleDb)) {
						warning("No alleleDb specified\n")
					}	else {
						if(is.null(alleleDb[[thisMarker]])) {
							warning(paste("No known alleles for",thisMarker))
							genotypeTable$allele.1 <- "noKnownAlleles"
							genotypeTable$allele.2 <- "noKnownAlleles"
							if(approxMatching) {	# need to match column names
								genotypeTable$allele.1.approx <- NA
								genotypeTable$allele.2.approx <- NA
							}
						} else  {
							if(is.null(resultObject@alleleDb[[thisMarker]])) {
								warning(paste("No variants for",thisMarker))
								genotypeTable$allele.1 <- NA
								genotypeTable$allele.2 <- NA
								if(approxMatching) {	# need to match column names
									genotypeTable$allele.1.approx <- NA
									genotypeTable$allele.2.approx <- NA
								}
							} else  {
								#varAlleleMap <- makeVarAlleleMap(alleleDb, resultObject@alleleDb, alleleMarkers=markerList, varMarkers=markerList)
								varAlleleMap <- makeVarAlleleMap(allele.variantMap=alleleDb[[thisMarker]], variant.variantMap=resultObject@alleleDb[[thisMarker]])
								genotypeTable$allele.1 <- varAlleleMap$knownAlleles[match(genotypeTable$varName.1, varAlleleMap$varNames)]
								genotypeTable$allele.2 <- varAlleleMap$knownAlleles[match(genotypeTable$varName.2, varAlleleMap$varNames)]

								if(approxMatching) {	# perform approx matching by BLAST
									cat("Attempting to find approximate matches using BLAST\n")
									varAlleleBlastMap <- makeVarAlleleBlastMap(allele.variantMap=alleleDb[[thisMarker]], variant.variantMap=resultObject@alleleDb[[thisMarker]])
									genotypeTable$allele.1.approx <- varAlleleBlastMap$subject[match(genotypeTable$varName.1, varAlleleBlastMap$query)]
									genotypeTable$allele.2.approx <- varAlleleBlastMap$subject[match(genotypeTable$varName.2, varAlleleBlastMap$query)]
								}
							}
						}

					}
				}
				callResults[[thisMarker]] <- new("genotypeCall", 
						projectName=resultObject@projectName, 
						runName=resultObject@runName,
						marker=thisMarker ,
						genotypeTable=genotypeTable,
						callMethod=method,
						callParameters=longParamList[[i]],
						mappedToAlleles=mapAlleles,
						alleleDbName="NeedToSetThis" )
			}
		}
	return(callResults)
}

### ver 0.14 have removed generics for callGenotypes to make mixed use of '...' and named parameters.
# generic method for callGenotypes.
## Make genotype calls
## 
## Apply a genotype call method to a table or list of tables of variant data such as the \code{markerSampleList} table of an \code{\link{mlgtResult}}.
##
## After \code{\link{mlgt}} has generated tables of the most common variants assigned in each marker/sample pair, an attempt can be made to call genotypes.
## This is kept separate because users might want to try different calling methods and have the option to map to a known set of alleles. Currently, only 
## one method is implemented (\emph{`custom'}).
## Methods:-
## \describe{
## 	\item{`callGenotypes.default'}{Three sequential steps for each marker/sample pair: 
## 		\enumerate{
##			\item {if the number of reads is less than \code{minTotalReads} the genotype is \emph{`tooFewReads'} }
##			\item {if the difference between the sum of counts of the top two variants and the count of the third most variant, expressed as proportion of total, is less than \code{minDiffToVarThree}, then the genotype is \emph{`complexVars'}}
##			\item {if the difference between the counts of top two variants, expressed as a proportion of the total, is greater than or equal to \code{minPropDiffHomHetThreshold}, then the genotype is \emph{HOMOZYGOTE}. Otherwise it is \emph{HETEROZYGOTE}. }  
##		}
## 	}
## }
## 
## 
## @param resultObject An object of class \code{\link{mlgtResult}}, as returned by \code{\link{mlgt}}
## @param table [Separate usage] A table of variant counts. 
## @param alleleDb A list of \code{\link{variantMap}} objects derived from known alleles. As made by 
## \code{\link{createKnownAlleleList}}
## @param method How to call genotypes. Currently only "callGenotypes.default" is implemented. Users can define
## their own methods as R functions (see the vignette).
## @param markerList For which of the markers do you want to call genotypes (default is all)?
## @param sampleList For which of the samples do you want to call genotypes (default is all)?
## @param mapAlleles FALSE/TRUE. Whether to map variants to db \option{alleleDb} of known alleles. 
##
## @return list of call results including the call parameters and a table of calls (class \code{\link{genotypeCall}}). If an mlgtResult object was supplied then a list of \code{\link{genotypeCall}} objects will be returned, each named by marker.
##
## @export
## @docType methods
## @aliases callGenotypes.mlgtResult
#setGeneric("callGenotypes",  function(resultObject, table, alleleDb=NULL, method="custom", minTotalReads=50, maxPropUniqueVars=0.8, 
#					minPropToCall=0.1, minDiffToVarThree=0.4,
#					minPropDiffHomHetThreshold=0.3, markerList=names(resultObject@markers),
#					sampleList=resultObject@samples, mapAlleles=FALSE, ...)
#	standardGeneric("callGenotypes")
#	) 


#setMethod("callGenotypes", signature(resultObject="missing", table="data.frame", alleleDb="ANY", method="ANY", minTotalReads="ANY", maxPropUniqueVars="ANY", 
#					minPropToCall="ANY", minDiffToVarThree="ANY",minPropDiffHomHetThreshold="ANY", markerList="ANY",
#					sampleList="ANY", mapAlleles="ANY"), 
#			definition=callGenotypes.default )


#setMethod("callGenotypes", signature(resultObject="mlgtResult", table="missing", alleleDb="ANY", method="ANY", minTotalReads="ANY", maxPropUniqueVars="ANY", 
#					minPropToCall="ANY", minDiffToVarThree="ANY",minPropDiffHomHetThreshold="ANY", markerList="ANY",
#					sampleList="ANY", mapAlleles="ANY"), 
#			definition=callGenotypes.mlgtResult )

#setMethod("callGenotypes", signature(resultObject="mlgtResult", table="missing"), 
#			definition=callGenotypes.mlgtResult )




# default for 'file' could be derived from project, run, marker attributes of genotypeCall.

writeGenotypeCallsToFile.genotypeCall <- function(genotypeCall, filePrefix="genoCalls", file=paste(filePrefix,genotypeCall@projectName,genotypeCall@runName,genotypeCall@marker,"tab",sep="."),
							writeParams=FALSE, appendValue=FALSE)  {
	if(writeParams)  {
		cat("# Genotype calls generated by R package 'mlgt' (", packageVersion("mlgt"),")", date(),"\n",file=file, append=appendValue)
		appendValue=TRUE
		cat("# Project:", genotypeCall@projectName,"\n", file=file, append=appendValue)
		cat("# Run", genotypeCall@runName,"\n", file=file, append=appendValue)
		cat("# Marker", genotypeCall@marker,"\n", file=file, append=appendValue)
		cat("# Call method:",  genotypeCall@callMethod,"\n", file=file, append=appendValue)		
		cat("# Call parameters:-\n", file=file, append=appendValue)
		for(i in 1:length(genotypeCall@callParameters)) {
			thisParam <- unlist(genotypeCall@callParameters[i])
			cat("#\t",names(thisParam),"=", thisParam,"\n", file=file, append=appendValue)
		}
		cat("# MappedToAlleles: ", genotypeCall@mappedToAlleles,"\n", file=file, append=appendValue)
		cat("# AlleleDb: ", genotypeCall@alleleDbName,"\n", file=file, append=appendValue)
	}	
	write.table(genotypeCall@genotypeTable, file=file, append=appendValue, row.names=F, quote=F, sep="\t")
	cat("Results written to",file,"\n")
	#return(invisible(1))
}


# Need function to write genotype calls as a single file.

writeGenotypeCallsToFile.list <- function(callList, filePrefix="genoCalls", file, singleFile=FALSE,writeParams=FALSE)  {
	if(singleFile)  {
		masterTable <- data.frame()
		for(i in 1:length(callList)) {
			masterTable <- rbind(masterTable , callList[[i]]@genotypeTable)
		}
		write.table(masterTable , file=file, row.names=F, quote=F, sep="\t")
		cat("Results written to",file,"\n")
	} else {
		#invisible(lapply(callList, FUN=function(x) writeGenotypeCallsToFile.genotypeCall(x,writeParams=writeParams)))
		cat(length(lapply(callList, FUN=function(x) writeGenotypeCallsToFile.genotypeCall(x,writeParams=writeParams, filePrefix=filePrefix))), "files written\n")
	}	

}

# DONE: allow write to single file when some markers mapped to alleles and some not. Consistent return of empty/NA/NaN columns?
# DONE: allow prefix to denote instance of filename - otherwise different genotype call sets will overwrite each other.
# generic method. Need to give defaults if defaults to be set in specific forms.
#' Write genotype calls to file
#' 
#' A genotype call table or a list of tables can be written to tab-delimited file(s). 
#'
#' This function is quite flexible and can output a single table of concatenated results or a series of individual files. 
#' Call parameters can be included above each table but be careful doing this when \code{singleFile=TRUE}
#' 
#' @param callList A \code{list} of genotypes calls.
#' @param genotypeCall Alternatively, supply a single table of genotype calls
#' @param filePrefix A prefix to add to the start of each file name. Useful to distinguish sets of genotype call results from same run.
#' @param file The file to write to. If none specified, function will attempt to make one. Ignored if \option{singleFile = TRUE}.
#' @param singleFile	FALSE/TRUE whether to concatenate results from a list of genotypeCalls
#' @param writeParams List call parameter values at top of file? Beware using this option when \option{singleFile = TRUE}
#' @param appendValue Used internally to concatenate results.
#'
#' @return Writes tables in the current working directory. 
#'
#' @rdname writeGenotypeCallsToFile-methods
#' @export
#' @aliases writeGenotypeCallsToFile.list,writeGenotypeCallsToFile.genotypeCall
#' @examples \dontrun{ 
#' data("mlgtResult", package="mlgt")
#' my.genoytpes <- callGenotypes(my.mlgt.Result)
#' writeGenotypeCallsToFile(my.genotypes)
#' }
setGeneric("writeGenotypeCallsToFile", function(callList, genotypeCall, filePrefix="genoCalls", file="", singleFile=FALSE, writeParams=FALSE, appendValue=FALSE)	standardGeneric("writeGenotypeCallsToFile")) 

#' @rdname writeGenotypeCallsToFile-methods
#' @aliases writeGenotypeCallsToFile,list,missing-method
setMethod("writeGenotypeCallsToFile", signature(callList="list", genotypeCall="missing", filePrefix="ANY", file="ANY", singleFile="ANY", 
									writeParams="ANY", appendValue="ANY"), 
				definition=writeGenotypeCallsToFile.list)

#' @rdname writeGenotypeCallsToFile-methods
#' @aliases writeGenotypeCallsToFile,missing,genotypeCall-method
setMethod("writeGenotypeCallsToFile", signature(callList="missing", genotypeCall="genotypeCall", filePrefix="ANY", file="ANY", singleFile="ANY", 
									writeParams="ANY", appendValue="ANY"), 
			definition=writeGenotypeCallsToFile.genotypeCall )

#writeGenotypeCallsToFile <- function(callList, genotypeCall, file="", singleFile=FALSE, writeParams=FALSE, appendValue=FALSE) attributes(callList)


########################## plotting
#' Plot BLAST statisitics for one marker
#' 
#' \code{\link{prepareMlgtRun}} produces several BLAST tables. It is instructive to plot the BLAST results and assess the performance
#'  of different markers. 
#'
#' This function is used to plot a series of histograms based on BLAST statistics.
#' 
#' @param blastTable The file of BLAST results.
#' @param subject The name of a single marker
#'
#' @return Plots three histograms based on the BLAST statistics 'Alignment length', 'Bit Score' and 'Percent Identity' 
#' @export
#' @seealso \code{\link{printBlastResultGraphs}}
inspectBlastResults <- function(blastTable, subject)  {

	#topHits <- getTopBlastHits(resultFile)
	#topHits$strand <- ifelse(topHits$s.end > topHits$s.start, 1,2)
	hitCount <- length(which(blastTable$subject == subject))
	if(hitCount > 0) { 
		#subject <- "DPA1_E2"
		breakValue <- max(10, 10^(floor(log10(hitCount))))	# favour a large number of breaks. At least 10.
		par(mfrow=c(1,3))
		hist(blastTable$ali.length[blastTable$subject == subject], breaks=breakValue, xlab="Alignment Length", main=subject)
		hist(blastTable$bit.score[blastTable$subject == subject], breaks=breakValue, xlab="Bit Score", main=subject)
		hist(blastTable$percent.id[blastTable$subject == subject], breaks=breakValue, xlab="% identity", main=subject)
	}	else  {
		warning(paste("No data for ", subject))
	}
}


#' Plot BLAST statistics for several markers to file
#' 
#' Plot the BLAST statistics easily for all markers of an \code{\link{mlgtResult}} object.
#'
#' 
#' 
#' @param designObject An object of class \code{\link{mlgtDesign}} which will contain the name of the blast results file \code{designObject@@markerBlastResults}
#' @param markerList Which markers to output. Defaults to \code{designObject@@markers}
#' @param fileName Defaults to "blastResultGraphs.pdf"
#'
#' @return Plots BLAST results to a pdf file.
#' @export
#' @seealso \code{\link{inspectBlastResults}}
printBlastResultGraphs <- function(designObject, markerList=designObject@markers, fileName="blastResultGraphs.pdf") {
	topHits <- getTopBlastHits(designObject@markerBlastResults)
	pdf(fileName,height=4)
	for(thisMarker in names(markerList))  {
		inspectBlastResults(topHits, thisMarker )
	}
	dev.off()
}



########## Multiple methods for plotGenotypeEvidence
# 
plotGenotypeEvidence.genotypeCall <- function(genotypeCall)  {
	genotypeTable <- genotypeCall@genotypeTable
	thisMarker	<- genotypeCall@marker
	# next three lines are a temp fix to retrieve default call parameters, which are not retained by callGenotypes.default(). May remove this function anyway.
	minTotalReads <- ifelse(exists("genotypeCall@callParameters[['minTotalReads']]"),genotypeCall@callParameters['minTotalReads'],50)
	minDiffToVarThree <- ifelse(exists("genotypeCall@callParameters[['minDiffToVarThree']]"),genotypeCall@callParameters['minDiffToVarThree'],0.4)
	minPropDiffHomHetThreshold <- ifelse(exists("genotypeCall@callParameters[['minPropDiffHomHetThreshold']]"),genotypeCall@callParameters['minPropDiffHomHetThreshold'],0.3)


	if(sum(genotypeTable$numbSeqs) < 1)  {
		warning(paste("No seqs to plot for",thisMarker), call.=F)
	} else {
		statusList <- as.factor(genotypeTable$status)
		pchList <- statusList
		levels(pchList) <- (1:nlevels(pchList ))
		#levels(pchList) <- 20+(1:nlevels(pchList ))


		par(mfrow=c(2,3))
		hist( genotypeTable$numbSeqs, breaks=20, main=thisMarker, xlab="numbSeqs"); abline(v=minTotalReads , lty=2)
		hist( genotypeTable$diffToVarThree, breaks=20, main=thisMarker, xlab="diffToVarThree", xlim=c(0,1)); abline(v=minDiffToVarThree , lty=2)
		hist(genotypeTable$propDiffHomHet, breaks=20, main=thisMarker, xlab="propDiffHomHet", xlim=c(0,1)) ; abline(v=minPropDiffHomHetThreshold , lty=2)

		plot(genotypeTable$diffToVarThree,genotypeTable$propDiffHomHet, main=thisMarker, xlab="diffToVarThree", ylab="propDiffHomHet",xlim=c(0,1), ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minDiffToVarThree , lty=2)
		legend("topleft", levels(as.factor(genotypeTable$status)), pch=as.numeric(levels(pchList)))
		plot(genotypeTable$numbSeqs,genotypeTable$diffToVarThree, main=thisMarker, xlab="numbSeqs", ylab="diffToVarThree", ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minDiffToVarThree , lty=2); abline(v=minTotalReads , lty=2)
		plot(genotypeTable$numbSeqs,genotypeTable$propDiffHomHet, main=thisMarker, xlab="numbSeqs", ylab="propDiffHomHet", ylim=c(0,1),pch=as.numeric(levels(pchList))[pchList]); abline(h=minPropDiffHomHetThreshold , lty=2); abline(v=minTotalReads , lty=2)
	}
}

plotGenotypeEvidence.genotypeCall.file <- function(genotypeCall, file)  {
	if(length(grep(".pdf$", file) ) < 1) {
		file <- paste(file,"pdf", sep=".")
	}
	pdf(file)
	plotGenotypeEvidence.genotypeCall(genotypeCall)
	dev.off()
	cat("Results output to", file, "\n")

}


# list must have a file specified for output
plotGenotypeEvidence.list <- function(callList, file) {
	if(length(grep(".pdf$", file) ) < 1) {
		file <- paste(file,"pdf", sep=".")
	}
	pdf(file)
	for(thisCall in callList) {
		#cat(thisCall@marker)
		plotGenotypeEvidence.genotypeCall(thisCall)
	}
	dev.off()
	cat("Results output to", file, "\n")	
	
}





#' Plot genotyping evidence
#' 
#' Plot the distributions of values used in calling genotypes. 
#'
#' Currently only makes sense with "custom" method. The resulting plots are
#' 	\enumerate{
#' 		\item {Histogram of the number of sequences assigned to each sample}
#'	 	\item {Histogram of diffToVarThree parameter. Used to decide whether to make the call}
#' 		\item {Histogram of propDiffHomHet parameter. Used to distinguish HOMOZYGOTES and HETEROZYGOTES}
#' 		\item {propDiffHomHet against diffToVarThree  }
#'	 	\item {diffToVarThree against number of sequences}
#' 		\item {propDiffHomHet against number of sequences}
#' }
#' 
#' @param callList A \code{list} of genotypes calls.
#' @param genotypeCall A single table of genotype calls
#' @param file The file to write to. 
#'
#' @return Creates six plots for each marker with a genotypeCall table. See \code{details}.
#'
#' @export
#' @seealso \code{\link{callGenotypes}}
#' @aliases plotGenotypeEvidence.list
#' @aliases plotGenotypeEvidence.genotypeCall.file
#' @aliases plotGenotypeEvidence.genotypeCall
#' @aliases plotGenotypeEvidence,missing,list,character-method
#' @aliases plotGenotypeEvidence,genoytpeCall,missing,missing-method
#' @aliases plotGenotypeEvidence,genoytpeCall,missing,character-method
#' @rdname plotGenotypeEvidence-methods
#' @examples \dontrun{ 
#' data("mlgtResult", package="mlgt")
#' my.genoytpes <- callGenotypes(my.mlgt.Result)
#' plotGenotypeEvidence(genotypeCall=my.genotypes[["DPA1_E2"]])
#' }
setGeneric("plotGenotypeEvidence", function(genotypeCall, callList, file) standardGeneric("plotGenotypeEvidence"))

#' @rdname plotGenotypeEvidence-methods
#' @aliases plotGenotypeEvidence,missing,list,character-method
setMethod("plotGenotypeEvidence", signature(genotypeCall="missing", callList="list", file="character"), definition=plotGenotypeEvidence.list)

#' @rdname plotGenotypeEvidence-methods
#' @aliases plotGenotypeEvidence,genotypeCall,missing,character-method
setMethod("plotGenotypeEvidence", signature(genotypeCall="genotypeCall", callList="missing", file="character"), definition=plotGenotypeEvidence.genotypeCall.file)

#' @rdname plotGenotypeEvidence-methods
#' @aliases plotGenotypeEvidence,genotypeCall,missing,missing-method
setMethod("plotGenotypeEvidence", signature(genotypeCall="genotypeCall", callList="missing", file="missing"), definition=plotGenotypeEvidence.genotypeCall)

#plotGenotypeEvidence <- function(callList, genotypeCall, file) attributes(genotypeCall)


#' Dump variants as fasta
#'
#' Output unique variants to one or more fasta files.
#'
#' This is a stop-gap function while I decide how best to handle output of full sequences. 
#' 
#' @param resultObject An object of class \code{\link{mlgtResult}} containing the sequence variants.
#' @param markers For which markers do you want to output sequences.
#' @param file An output file name. If not supplied, one is created.
#' @param singleFile Whether to output results for all markers to a single file or to one file per marker.
#'
#' @return Writes fasta files in the current directory. 
#' @export
dumpVariantMap.mlgtResult <- function(resultObject, markers=names(resultObject@markers),
			 file=paste(resultObject@projectName,resultObject@runName,"seqDump",sep="."),
			singleFile=TRUE) {

	file <- sub("\\.fasta","",file)
	writeOption <- ifelse(singleFile,"a","w")	# used by write.fasta. "a" = append
	baseFileName <- file	
	for(thisMarker in markers)  {
		useFile <- ifelse(singleFile,paste(baseFileName,"fasta",sep="."),paste(baseFileName,thisMarker,"fasta",sep="."))	
		theSeqs <- names(resultObject@alleleDb[[thisMarker]]@variantMap)
		theNames <-  as.character(resultObject@alleleDb[[thisMarker]]@variantMap)
		write.fasta(lapply(theSeqs,s2c), theNames, file.out=useFile , open=writeOption )
	}
}


#dumpVariantMap.variantMap


## TODO docs, internal
stripGapColumns <- function(alignment)  {
	gap_index <- which(con(alignment, method="threshold",threshold=(1-1e-07)) == '-')		# bug in seqinr: threshold =1 returns NA for all.
	if(length(gap_index) < 1)  {	# nothing to strip
		if(exists("verbose")) cat("No gap columns to remove ")
		return(alignment)
	}
	alignMatrix <- as.matrix(alignment)[,-c(gap_index)]
	if(exists("verbose")) cat(paste("removed",length(gap_index),"gap columns "))
	deGapSeqs <- apply(alignMatrix,1,c2s)
	deGapAlign <- as.alignment(alignment$nb,alignment$nam,deGapSeqs)
	return(deGapAlign)
}


## TODO docs, internal
## supply an UN-PACKKED alignment (including duplicated sequences)
errorCorrect.alignment <- function(alignment, correctThreshold=0.01)  {
	#alignment.corrected <- alignment
	minDepth <- ceiling(1/correctThreshold)	# no point attempting if less depth than this.
	if(alignment$nb < minDepth) {
		warning(paste("No correction possible with depth of", alignment$nb, "and correction threshold of", correctThreshold, "\n"))
		return(alignment)
	}
		thisProfile <- consensus(alignment, method="profile")
		thisConsensus <- con(alignment, method="threshold", threshold=correctThreshold )
		totalSeqs <- alignment$nb
		totalLetters <- nrow(thisProfile)
		mafList <- apply(thisProfile,  2, FUN=function(x) (sort(x)[totalLetters-1] / totalSeqs)) 
		correct_index <- intersect(which(mafList < correctThreshold), which(mafList > 0))
		#remove NAs from index.
		correct_index <- correct_index[!is.na(thisConsensus[correct_index])]
	
		
		alignment.matrix <- as.matrix.alignment(alignment)
		alignment.matrix[,c(correct_index)] <- t(apply(alignment.matrix,1,FUN=function(x)  x[c(correct_index)]  <- thisConsensus[correct_index]))
		seqList.correct <- apply(alignment.matrix,1,c2s)
		alignment.corrected <- as.alignment(nb=length(seqList.correct),nam=names(seqList.correct),seq=seqList.correct)
		alignment.corrected <- stripGapColumns(alignment.corrected)
		return(alignment.corrected)
}



## TODO docs, incorporate into mlgt() lots of duplicated code, README
##
## wrapper function to perform errorCorrect on an existing mlgtResult object
## this is basically mlgt() without any of the assignment and alignment and a call to errorCorrect()
## Examples:-
## my.corrected.mlgt.Result <- errorCorrect.mlgtResult(my.mlgt.Result)
## my.corrected.mlgt.Result <- errorCorrect.mlgtResult(my.mlgt.Result,correctThreshold=0.05)
## alignReport(my.mlgt.Result, fileName="alignReportOut.pdf", method="profile")
## alignReport(my.corrected.mlgt.Result, fileName="alignReportOut.corrected.pdf", method="profile")
## alignReport(my.mlgt.Result, fileName="alignReportOut.hist.pdf", method="hist")
## alignReport(my.corrected.mlgt.Result, fileName="alignReportOut.hist.corrected.pdf")
## my.genotypes <- callGenotypes(my.mlgt.Result)
## my.corrected.genotypes <- callGenotypes(my.corrected.mlgt.Result)
## my.genotypes[[thisMarker]]@genotypeTable
## my.corrected.genotypes[[thisMarker]]@genotypeTable
## error correction has little effect on the the sample dataset even with high threshold of 0.05. Most samples with > 100 seqs are called correctly anyway.
errorCorrect.mlgtResult  <- function(mlgtResultObject, correctThreshold=0.01)  {



	markerSampleList <- list()
	runSummaryTable <- data.frame()
	alleleDb <- list()
	varCountTableList <- list()


	###ITERATIONS
	for(thisMarker in names(mlgtResultObject@markers)) {


	cat(paste(thisMarker,"\n"))
	#thisMarker <- "DQA1_E2"

	## might need to combine all these to return a single item.
	summaryList <- list()
	summaryTable <- data.frame()
	markerSequenceCount <- list("noSeq"=0)		#  BUG? requires some data otherwise won't sum properly with localSequenceCount.
	alleleList <- list() 
	variantList <- list()
	alleleCount <- 1
	markerSeq <- unlist(getSequence(mlgtResultObject@markers[[thisMarker]],as.string=T))
	varCountTableList[[thisMarker]] <- data.frame()

	thisTable <- mlgtResultObject@varCountTables[[thisMarker]]
	
	for(thisSample in mlgtResultObject@samples) {
			cat(paste(thisSample ," "))
		## need to keep data from mlgtResultObject where no correction is made, but create new data where a correction is made.

			#testPairSeqList <- intersect(pairedSampleMap[[thisSample]],union(fMarkerMap[[thisMarker]], rMarkerMap[[thisMarker]]))



		seqTable <- data.frame()
		localAlleleNames <- c("NA","NA","NA")
		localAlleleFreqs <- c(0,0,0)

		## go through all seq's mapped to this marker/sample pair.
		## extract the corresponding sequence delimited by the top blast hits on the primers.  IS THIS THE BEST WAY?
		##		Simple improvement: minimum blast hit length to primer to keep. 

		## internal Function

		recordNoSeqs <- function(summaryTable)  {		# to record no seqs before skipping out. 
				summaryRow <- data.frame(marker=thisMarker, sample=thisSample, numbSeqs=0,numbVars=0,
					varName.1="NA", varFreq.1= 0,
					varName.2="NA", varFreq.2= 0,
					varName.3="NA", varFreq.3= 0)
				summaryTable <- rbind(summaryTable, summaryRow)
				return(summaryTable)
		}

			if(is.na(match(thisSample,names(thisTable)))) {
				summaryTable  <- recordNoSeqs(summaryTable)
				next;
			}

			valueIndex <- !is.na(thisTable[,thisSample])
			seqCounts <- thisTable[valueIndex,thisSample]
			## important to 'unpack' the alignment so that each sequence occurs the correct number of times.
			sampleSeqs <- rep(row.names(thisTable)[valueIndex ], seqCounts)
			thisAlign <- as.alignment(sum(seqCounts), sampleSeqs, sampleSeqs)
			cat(paste(length(unique(thisAlign$seq)),"/", thisAlign$nb,"unique seqs in original alignment, "))

			newAlign <- errorCorrect.alignment(thisAlign, correctThreshold)
			
			cat(paste(length(unique(newAlign$seq)),"/", newAlign$nb,"unique seqs in new alignment, "))
			## DONE: repack the corrected alignment re-attribute the allele names. Update the markerSampleTable
			varCountTable <-  as.data.frame(table(unlist(newAlign$seq)),stringsAsFactors=FALSE)

			seqTable <- data.frame(alignedVar=varCountTable$Var1, count=as.numeric(varCountTable$Freq),stringsAsFactors=FALSE)
			seqTable$var <- gsub("-","",seqTable$alignedVar)
			seqTable <- seqTable[order(seqTable$count,decreasing=T),]
		#dim(seqTable)
		# if no sequences returned, nothing to process. 
		if(nrow(seqTable) < 1 )  {
			summaryTable  <- recordNoSeqs(summaryTable)
			cat(" no variants\n")
			#summaryList[[thisMarker]][[thisSample]] <- NA	
			next ;		# go to next sample.
		}

		# store sequence and count of sequence as alignedVar (needed for alignment report)
		varCountTableList[[thisMarker]][seqTable$alignedVar,thisSample] <- seqTable$count

		## test if variants are novel. 
		## Give allele names?  
		## Do with first three for now. 

		alToRecord <- min(3,nrow(seqTable))
		if(alToRecord > 0)  {
			for (a in 1:alToRecord )  {
				if(is.null(variantList[[seqTable$var[a]]]))  {    	# novel
					alleleName <- paste(thisMarker, alleleCount,sep=".")	
					variantList[[seqTable$var[a]]] <- alleleName
					localAlleleNames[a] <- alleleName 
					localAlleleFreqs[a] <- seqTable$count[a]
					alleleCount <- alleleCount + 1
				} else  {										# pre-existing alllele
					localAlleleNames[a] <- variantList[[seqTable$var[a]]]
					localAlleleFreqs[a] <- seqTable$count[a]		
				}
			}
		}


		# compile stats

		if(nrow(seqTable) >0 )  {	# cannot allow assignment from empty list as messes up class of list for remaining iterations
			summaryList[[thisMarker]] <- list()
			summaryList[[thisMarker]][[thisSample]] <- seqTable
		}

		summaryRow <- data.frame(marker=thisMarker, sample=thisSample, numbSeqs=sum(seqTable$count),numbVars=nrow(seqTable),
					varName.1=localAlleleNames[1], varFreq.1= localAlleleFreqs[1],
					varName.2=localAlleleNames[2], varFreq.2= localAlleleFreqs[2],
					varName.3=localAlleleNames[3], varFreq.3= localAlleleFreqs[3])
		summaryTable <- rbind(summaryTable, summaryRow)

		#sequence count across samples? 
		# need to sum from summaryTable or from summaryList.
		#markerSequenceCount <- 
		#as.list(colSums(merge(m, n, all = TRUE), na.rm = TRUE))  # not working
		localSequenceCount <- as.list(seqTable$count)
		names(localSequenceCount) <- seqTable$var
		markerSequenceCount   <- as.list(colSums(merge(markerSequenceCount  , localSequenceCount,  all = TRUE), na.rm = TRUE))
		# might need to instantiate the markerSequenceCount if empty. 

			cat("\n")

	}  # end of sample loop

	markerSampleList[[thisMarker]] <- summaryTable
	## DONE: min(nchar(names(variantList))) throws warning when no variants in list. (Inf/-Inf)
	minVarLength <- ifelse(length(markerSequenceCount) < 1, NA, min(nchar(names(markerSequenceCount))) )
	maxVarLength <- ifelse(length(markerSequenceCount) < 1, NA, max(nchar(names(markerSequenceCount))) )
	minAleLength <- ifelse(length(variantList) < 1, NA, min(nchar(names(variantList))) )
	maxAleLength <- ifelse(length(variantList) < 1, NA, max(nchar(names(variantList))) )

	runSummaryRow <- data.frame(marker=thisMarker, assignedSeqs=sum(summaryTable$numbSeqs), assignedVariants=sum(summaryTable$numbVars), 
					minVariantLength=minVarLength, 
					maxVariantLength=maxVarLength,
					minAlleleLength=minAleLength, maxAlleleLength=maxAleLength )
	runSummaryTable <- rbind(runSummaryTable, runSummaryRow)
	if(length(variantList) > 0)  {
		# This line replaced. Not entirely tested the repurcussions. e.g. makeVarAlleleMap()?
		#alleleDb[[thisMarker]] <- variantList   # LATER: separate lists for alleles and variants? 
		#alleleDb[[thisMarker]] <- list(reference=as.SeqFastadna(markerSeq, thisMarker), alleleMap=variantList, inputAlleleCount = length(unlist(variantList)), uniqueSubAlleleCount=length(variantList))
		alleleDb[[thisMarker]] <- new("variantMap", reference=as.SeqFastadna(markerSeq, thisMarker), 
							variantSource=paste(mlgtResultObject@projectName, mlgtResultObject@runName,sep="."),
							variantMap=variantList, inputVariantCount = length(unlist(variantList)), uniqueSubVariantCount=length(variantList))
	}

}  # end of marker loop
	
	localMlgtResult <- new("mlgtResult", mlgtResultObject,  runSummaryTable=runSummaryTable , alleleDb=alleleDb, markerSampleList=markerSampleList,
						varCountTables=varCountTableList)
	return(localMlgtResult)

} 


#' Alignment error correction
#' 
#' Correct very low frequency site variants.
#' 
#' You may want to alter some of the sequences if you believe that sequences at very low frequency 
#' (within the set of sequences from a marker/sample pair) represent sequencing errors.
#' \code{errorCorrect()} is implemented as an additional step after running \code{\link{mlgt}}, however, it is recommended to include error correction within \code{\link{mlgt}} using the errorCorrect=TRUE option.
#' Using \code{\link{alignReport}} beforehand may help you decide whether to do this. 
#' 
#' @param mlgtResultObject An object of class \code{\link{mlgtResult}}
#' @param correctThreshold The maximimum Minor Allele Frequency (MAF) at which variants will be corrected. 
#' 
#' @return A new \code{\link{mlgtResult}} object with errors 'corrected'
#' @rdname errorCorrect-methods
#' @export 
#' @seealso \code{\link{alignReport}}
#'
setGeneric("errorCorrect", function(mlgtResultObject, alignment, correctThreshold=0.01) standardGeneric("errorCorrect")) 

#' @rdname errorCorrect-methods
#' @aliases errorCorrect,missing,list-method
setMethod("errorCorrect", signature(mlgtResultObject="missing",alignment="list", correctThreshold="ANY"), 
				definition=errorCorrect.alignment)

#' @rdname errorCorrect-methods
#' @aliases errorCorrect,mlgtResult,missing-method
setMethod("errorCorrect", signature(mlgtResultObject="mlgtResult", alignment="missing", correctThreshold="ANY"), 
				definition=errorCorrect.mlgtResult)



## Generate stats per site along the alignments. WITHIN a marker/sample pair.
## DONE: Docs needed, add e.g. to README
## This function is a bit weird in that it collects a table of info and, optionally, generates some graphs.
## What if people want to export the tables of results to file AND the images?
## EXAMPLES:- 
## alignReport(my.mlgt.Result,markers=thisMarker, samples=thisSample, method="profile")
## alignReport(my.mlgt.Result,markers=thisMarker, samples=thisSample, method="hist", fileName="testOutHist")
## alignReport(my.mlgt.Result, method="profile", fileName="testOutMultiProfile")
## alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-17", method="profile")
## alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-1", method="hist")
## alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-22", method="hist") # good example where change would be useful
## alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-1", method="profile", correctThreshold=0.02)
## my.alignReport <- alignReport(my.mlgt.Result)
#' Report on alignment
#' 
#' Inspect site frequency spectra for alignments.
#' 
#' Produce different kinds of reports to assess quality of data for each marker/sample pair. 
#' Can be a good way to assess whether \code{\link{errorCorrect}} should be applied.
#' 
#' 
#' @param mlgtResultObject an object of class \code{\link{mlgtResult}}
#' @param markers Which markers to output
#' @param samples Which samples to output
#' @param correctThreshold A hypothetical level at which you migth correct low frequence variants. Default = 0.01.
#' @param consThreshold (1- correctThreshold)
#' @param profPlotWidth How many residues to plot in \code{profile} mode. Default=60.
#' @param fileName Give a filename to export result to (pdf).
#' @param method One of c("table", "profile", "hist"). "hist" plot a histogram of MAF frequencies. "profile" plots a coloured barplot represnting the allele frequencies at each site.
#' @param warn Issue warnings (default = TRUE)
#' 
#' @return A data frame for each marker listing site statistics.
#' @seealso \code{\link{errorCorrect}}
#' @export 
#' @examples \dontrun{ 
#' data("mlgtResult", package="mlgt")
#' alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-22", method="profile")
#' alignReport(my.mlgt.Result,markers="DPA1_E2", samples="MID-22", method="hist")
#' }
alignReport <- function(mlgtResultObject, markers=names(mlgtResultObject@markers), samples=mlgtResultObject@samples,
		correctThreshold = 0.01,  consThreshold = (1 - correctThreshold), profPlotWidth = 60, fileName=NULL, method="table", warn=TRUE)  {

	# need method for both plots (save processing time) but how to do without generating profiles twice.
	# need to tidy and label profile plots.

	reportList <- list()
	if(!is.null(fileName)) {
		if(length(grep(".pdf$", fileName))  < 1 & !is.null(fileName)) {
				fileName <- paste(fileName,"pdf", sep=".")
			}
		pdf(fileName)
	}

	colTable <- data.frame(profChar = c("-","a","c","g","t"),
					profCol = c("grey","green", "blue", "yellow", "red"))

	for(thisMarker in markers)  {
		thisTable <- mlgtResultObject@varCountTables[[thisMarker]]
		cat(paste(thisMarker,"\n"))
		if(nrow(thisTable) < 1)  {		# drops this marker if no data. Might be better to return and empty table?
			if(warn) { 
					warning(paste("No data for",  thisMarker)) 
				}
			#reportTable[thisSample,c("numbSeqs","numbVars")] <- 0
			#reportTable[thisSample,c("alignLength","invar.sites","mafBelowThreshold","mafAboveThreshold")] <- NA
			next;
		}
		reportTable <- data.frame()
		for(thisSample in samples) {
			if(is.na(match(thisSample,names(thisTable)))) {
				if(warn) { 
					warning(paste("No variants for", thisSample, "and", thisMarker)) 
				}
				#reportTable[thisSample,c("invar.sites","mafBelowThreshold","mafAboveThreshold")] <- NA
				reportTable[thisSample,"numbSeqs"] <- 0
				reportTable[thisSample,"numbVars"] <- 0
				reportTable[thisSample,"alignLength"] <- NA
				reportTable[thisSample,"invar.sites"] <- NA
				reportTable[thisSample,"mafAboveThreshold"] <- NA
				reportTable[thisSample,"mafBelowThreshold"] <- NA
				next;
			}
			#cat(paste(thisSample ," "))
			valueIndex <- !is.na(thisTable[,thisSample])
			seqCounts <- thisTable[valueIndex,thisSample]
			sampleSeqs <- rep(row.names(thisTable)[valueIndex ], seqCounts)
			thisAlign <- as.alignment(sum(seqCounts), sampleSeqs, sampleSeqs)

			if(thisAlign$nb < 2)  {	# too few sequences to plot.
				#reportTable[thisSample,c("invar.sites","mafBelowThreshold","mafAboveThreshold")] <- NA
				reportTable[thisSample,"numbSeqs"] <- 1
				reportTable[thisSample,"numbVars"] <- 1
				reportTable[thisSample,"alignLength"] <- nchar(thisAlign$seq[1])
				reportTable[thisSample,"invar.sites"] <- NA
				reportTable[thisSample,"mafAboveThreshold"] <- NA
				reportTable[thisSample,"mafBelowThreshold"] <- NA
				next; 
			}
			thisProfile <- consensus(thisAlign , method="profile")
			thisConsensus <- con(thisAlign, method="threshold", threshold=consThreshold)
			totalSeqs <- sum(seqCounts)
			totalLetters <- nrow(thisProfile)
			mafList <- apply(thisProfile,  2, FUN=function(x) (sort(x)[totalLetters-1] / totalSeqs)) 

			reportTable[thisSample, "numbSeqs"] <- totalSeqs 
			reportTable[thisSample, "numbVars"] <- length(seqCounts)
			reportTable[thisSample, "alignLength"] <- ncol(thisProfile)			
			reportTable[thisSample, "invar.sites"] <- sum(mafList == 0)
			## variable sites with minor allele > correction threshold.
			reportTable[thisSample, "mafAboveThreshold"] <- sum(mafList >= correctThreshold )
			## varaible sites with minor alleles < correction threshold.
			reportTable[thisSample, "mafBelowThreshold"] <- reportTable[thisSample, "alignLength"] - (reportTable[thisSample, "invar.sites"] + reportTable[thisSample, "mafAboveThreshold"])
	
			if(method=="profile")  {

				profColours <-  as.character(colTable$profCol[match(row.names(thisProfile),colTable$profChar)])
				## splits the plotting across so many lines. Buffers final plot to constant length.
				profLen <- ncol(thisProfile)
				n_plot <- ceiling(profLen / profPlotWidth )
				plotProfile <- thisProfile
				plotConsensus <-  toupper(thisConsensus) 
				remainder <- profLen %% profPlotWidth 

				if(remainder > 0) {
					extLen <- profLen + (profPlotWidth - remainder)
					profExtension <- matrix(0,nrow=nrow(thisProfile), ncol=(profPlotWidth - remainder), 
								dimnames=list(row.names(thisProfile),c((profLen+1): extLen)))
					plotProfile <- cbind(thisProfile,profExtension)
					plotConsensus <- c(toupper(thisConsensus), rep("",remainder)) 
				}

				old.o <- par(mfrow=c((n_plot+1),1), mar=c(2,4,2,2))
				plot.new()	
				title(paste(thisMarker, thisSample, sep=" : "), line=0)		
				title("", line=0,adj=0,	sub=paste(c("Total sites","Invariant sites","MAF above threshold","MAF below threshold"),reportTable[thisSample,3:6],sep=": ",collapse="\n") )
				
				legend("right", legend=toupper(colTable$profChar),fill=as.character(colTable$profCol), horiz=T)
				for(i in 1:n_plot) {
					start <- ((i-1)*profPlotWidth) + 1
					#index <- start:(min(start+profPlotWidth, profLen))

					index <- start:((start+profPlotWidth)-1)
					barplot(plotProfile[,index ], col=profColours , names.arg=toupper(plotConsensus[index]) )
				}
				par(old.o)
			}

			if(method=="hist")  {
				#if(!is.null(fileName)) pdf(fileName)
				if(sum(mafList > 0) > 0) {
					hist(mafList[mafList > 0], breaks=200, xlim=c(0,0.5), xlab="Site-specific minor allele frequency", sub="non-zero values only",main=paste(thisMarker, thisSample, sep=":")) 
					abline(v=correctThreshold, lty=2)
				}	
			}
			
		}			
		reportList[[thisMarker]] <- reportTable
	}
	if(!is.null(fileName)) {
		dev.off()
		cat(paste("Alignment figures(s) plotted to", fileName,"\n"))
	}
	#print(reportList)
	return(reportList)
}



## DONE docs, README entry
## Examples
## dumpVariants(my.mlgt.Result,fileSuffix="variantDump.fasta")
## dumpVariants(my.corrected.mlgt.Result, markers="FLA_DRB", samples="cat348" )
## dumpVariants(my.corrected.mlgt.Result,fileSuffix="corrected.0.05.variantDump.fasta")
## dumpVariants(my.corrected.mlgt.Result,fileSuffix="corrected.0.05.variantDump.unique.fasta", uniqueOnly=T)
#' Print sequence to file
#' 
#' A function to output all sequences or just unique sequences to a fasta file
#' 
#' The sequence variants stored within an object of class \code{\link{mlgtResult}}
#' are not very easy to extract. This function will output all variants or all 
#' variant for specific markers and samples into fasta files. Users can select to only
#' output unique sequences or the full alignment including duplicated sequences.
#' One file will be created for each marker/sample pair. 
#' 
#' @param mlgtResultObject an object of class \code{\link{mlgtResult}}
#' @param markers Which markers to output
#' @param samples Which samples to output
#' @param fileSuffix Add a common suffix to the file names. Usefull for keeping track of different sets of sequences.
#' @param uniqueOnly Only output single copy of each sequence. A count for each sequence are appended to the names.
#' 
#' @export
#' 
dumpVariants <- function(mlgtResultObject, markers=names(mlgtResultObject@markers),
			 samples=mlgtResultObject@samples, fileSuffix="variantDump.fasta", uniqueOnly=FALSE) {
	for(thisMarker in markers)  {
		thisTable <- mlgtResultObject@varCountTables[[thisMarker]]
		for(thisSample in samples)  {
			
			if(is.na(match(thisSample,names(thisTable)))) {
				warning(paste("Nothing to output for", thisMarker, thisSample))
				next;
			}
			fileName <- paste(thisMarker,thisSample,fileSuffix, sep=".")
			valueIndex <- !is.na(thisTable[,thisSample])
			seqCounts <- thisTable[valueIndex,thisSample]
			if(uniqueOnly) {		# export one copy of each sequence, append count to name line
				sampleSeqs <- row.names(thisTable)[valueIndex]
				seqNames <- paste(sampleSeqs,seqCounts)
				thisAlign <- as.alignment(length(seqCounts), seqNames, sampleSeqs)
			} else {			# export copies of seqs as per counts.
				## important to 'unpack' the alignment so that each sequence occurs the correct number of times.
				sampleSeqs <- rep(row.names(thisTable)[valueIndex ], seqCounts)
				thisAlign <- as.alignment(sum(seqCounts), sampleSeqs, sampleSeqs)
			}
			write.fasta(sequences=lapply(thisAlign$seq,s2c), names=thisAlign$nam, file.out=fileName)
			cat(paste("Variants written to", fileName,"\n"))
		}
	}
}



## used by combineMlgtResults()
mergeMlgtResults.complex <- function(result1, result2) {
	master <- result1
	# need to check if samples shared between results.
	# If not, then perform complex join. Combining results where approriate.
	# If they are, then test if marker/sample pairs shared. If yes, then stop(), otherwise perform complex join.
	master@markers <- c(master@markers,result2@markers)
	master@markers <- master@markers[!duplicated(master@markers)]
	newRunSummaryTable <- data.frame()
	for(thisMarker in names(master@markers)) {
		markerSampleOverlap <- intersect(master@markerSampleList[[thisMarker]]$sample,result2@markerSampleList[[thisMarker]]$sample)
		if(length(markerSampleOverlap) > 0)  {
			# STOP!
			stop(paste("Cannot join results sharing marker results for same samples\n",thisMarker,markerSampleOverlap,"\n"))
		}

		# need to combine alleleDBs and variantMaps first so that new consistent allele names can propagate to the markerSampleList
		#master@alleleDb <- as.list(merge(master@alleleDb, result2@alleleDb))
		#master@alleleDb <- mergeAlleleDbs(master@alleleDb)
		masterAlleleTable <- data.frame(seq=names(unlist(master@alleleDb[[thisMarker]]@variantMap)),
								masterName=as.character(unlist(master@alleleDb[[thisMarker]]@variantMap)), stringsAsFactors=F )
		alleleTable2 <- data.frame(seq=names(unlist(result2@alleleDb[[thisMarker]]@variantMap)),
								alleleName=as.character(unlist(result2@alleleDb[[thisMarker]]@variantMap)), stringsAsFactors=F )
		masterMatchTable <- merge(masterAlleleTable, alleleTable2, by="seq", all=T)
		# some alleles will be new and these need to be added to master table and given new allele names
		maxAllele <- max(na.omit(as.numeric(sub(paste(thisMarker,".",sep=""),"",masterMatchTable$masterName))))
		newAlleleCount <-sum(is.na(masterMatchTable$masterName))
		masterMatchTable$masterName[is.na(masterMatchTable$masterName)] <- paste(thisMarker, (1:newAlleleCount + maxAllele), sep=".")
	
		#create new variantMap from masterMatchTable
		newVarMap <- as.list(masterMatchTable$masterName)
		names(newVarMap) <- masterMatchTable$seq
		master@alleleDb[[thisMarker]]@variantMap <- newVarMap
	
		# update allele Names in result2 table
		nameCols <- grep("varName", names(result2@markerSampleList[[thisMarker]]))
		for(n in nameCols) {
			result2@markerSampleList[[thisMarker]][,n] <- masterMatchTable$masterName[match(result2@markerSampleList[[thisMarker]][,n] , masterMatchTable$alleleName)]
		}
		master@markerSampleList[[thisMarker]] <- rbind(master@markerSampleList[[thisMarker]],result2@markerSampleList[[thisMarker]])
	
		# Combine result across varCountTables.
		#master@varCountTables <- as.list(merge(master@varCountTables, result2@varCountTables))
		# varCountTables[[thisMarker]][seqTable$alignedVar,thisSample] <- seqTable$count
		frame2 <- result2@varCountTables[[thisMarker]]
		for(thisCol in names(frame2)) {
			index <- which(!is.na(frame2[,thisCol] ))
			counts <- frame2[index,thisCol]
			names <- row.names(frame2)[index]
			master@varCountTables[[thisMarker]][names,thisCol] <- counts
		}

		index1 <- match(thisMarker, master@runSummaryTable$marker)
		index2 <- match(thisMarker, result2@runSummaryTable$marker)

		runSummaryRow <- data.frame(marker=thisMarker, assignedSeqs=sum(master@markerSampleList[[thisMarker]]$numbSeqs), assignedVariants=sum(master@markerSampleList[[thisMarker]]$numbVars), 
					minVariantLength=min(master@runSummaryTable$minVariantLength[index1], result2@runSummaryTable$minVariantLength[index2] ), 
					maxVariantLength=max(master@runSummaryTable$maxVariantLength[index1], result2@runSummaryTable$maxVariantLength[index2] ),
					minAlleleLength=min(master@runSummaryTable$minAlleleLength[index1], result2@runSummaryTable$minAlleleLength[index2] ), 
					maxAlleleLength=max(master@runSummaryTable$maxAlleleLength[index1], result2@runSummaryTable$maxAlleleLength[index2] ) )
		newRunSummaryTable <- rbind(newRunSummaryTable,runSummaryRow )
	}	
	master@samples <- union(master@samples,result2@samples)
	master@runSummaryTable <- newRunSummaryTable 
	# keep the following unless different between results.
	if(!identical(master@fTags,result2@fTags)) {	master@fTags <- list() }	# @fTags. List of class SeqFastadna
	if(!identical(master@rTags,result2@rTags)) {	master@rTags<- list() }		# @rTags. List of class SeqFastadna
	if(!identical(master@inputFastaFile,result2@inputFastaFile)) {	master@inputFastaFile <- '' }		# @inputFastaFile    
	if(!identical(master@markerBlastResults,result2@markerBlastResults)) {	master@markerBlastResults <- '' }		# @markerBlastResults
	if(!identical(master@fTagBlastResults,result2@fTagBlastResults)) {	master@fTagBlastResults <- '' }		# @fTagBlastResults  
	if(!identical(master@rTagBlastResults,result2@rTagBlastResults)) {	master@rTagBlastResults <- '' }		# @rTagBlastResults  
	if(!identical(master@runName,result2@runName)) {	master@runName<- 'CombinedResults' }			# @runName
	# @projectName is taken as the first result. May need to change this.
	return(master)
}

## used by combineMlgtResults()
mergeMlgtResults.simple <- function(result1, result2) {
	master <- result1
	master@markers <- c(master@markers,result2@markers)
	master@markers <- master@markers[!duplicated(master@markers)]
	master@markerSampleList <- c(master@markerSampleList, result2@markerSampleList)
	master@markerSampleList <- master@markerSampleList[!duplicated(master@markerSampleList)]
	master@alleleDb <- c(master@alleleDb, result2@alleleDb)
	master@alleleDb <- master@alleleDb[!duplicated(master@alleleDb)]
	master@varCountTables <- c(master@varCountTables, result2@varCountTables)
	master@varCountTables <- master@varCountTables[!duplicated(master@varCountTables)]
	master@samples <- union(master@samples,result2@samples)
	master@runSummaryTable <- rbind(master@runSummaryTable, result2@runSummaryTable)
	# keep the following unless different between results.
	if(!identical(master@fTags,result2@fTags)) {	master@fTags <- list() }	# @fTags. List of class SeqFastadna
	if(!identical(master@rTags,result2@rTags)) {	master@rTags<- list() }		# @rTags. List of class SeqFastadna
	if(!identical(master@inputFastaFile,result2@inputFastaFile)) {	master@inputFastaFile <- '' }		# @inputFastaFile    
	if(!identical(master@markerBlastResults,result2@markerBlastResults)) {	master@markerBlastResults <- '' }		# @markerBlastResults
	if(!identical(master@fTagBlastResults,result2@fTagBlastResults)) {	master@fTagBlastResults <- '' }		# @fTagBlastResults  
	if(!identical(master@rTagBlastResults,result2@rTagBlastResults)) {	master@rTagBlastResults <- '' }		# @rTagBlastResults  
	if(!identical(master@runName,result2@runName)) {	master@runName<- 'CombinedResults' }			# @runName
	# @projectName is taken as the first result. May need to change this.
	return(master)
}

## I imagine most of the time, combining results can be done after mlgt is run by simple concatenation of genotype tables. 
## This might be the case if certain markers are rerun. 
## However, there are instances where combination within mlgt is desirable. e.g. after a parallel run. 
## Re-runs of certain marker samples could also be accomodated, but care would be have to be taken to update the correct results.
## 1. to combine and summarise results across runs.
## 2. to re-combine results after running parallel processing.
## Lots of the merging should be in common.
## Need checks that samples/markers are not duplicated (option to rename?)
## Needs to be fast enough for point 2 above to be worthwhile.
## Even if markers don't overlap, need to check that samples and mids are the same or different.
## Or do I? MIDs should be allowed to differ between runs. 
#' Combine two or more mlgtResult objects
#' 
#' Combine results from one or more runs, or combine partial results after a parallel job.
#' 
#' In some cases, you may want to combine multiple \code{\link{mlgtResult}} objects
#' into a single object. 
#' Can combine results using the same markers as long as the samples used have different names between results.
#' Can combine results using different sets (subsets) of markers.
#' Will fail if the same marker/sample combination appears in more than one \code{\link{mlgtResult}}.
#' Can be used to recombine the list of result obtained by running \code{\link{mlgt}} in parallel
#' on subsets of the full marker list.
#'
#'
#' @param resultList A list of objects of class \code{\link{mlgtResult}}
#' @param projectName Do you want to provide your own projectName
#' @param runName Do you want to provide your own runName
#' 
#' @return An object of class \code{\link{mlgtResult}}
#' @export
#'
combineMlgtResults <- function(resultList,projectName=resultList[[1]]@projectName, runName="combinedMlgtResults")  {
	# set first result as master
	master <- resultList[[1]]

	for(i in 2:length(resultList)) {
		# cycle through remaining results adding each to master
		# check that sample/marker combinations do not overlap
		#Any overlap in markers?
		markerOverlap <- intersect(names(master@markers),names(resultList[[i]]@markers))
			
		if(length(markerOverlap) > 0) {	# overlapping markers, test for sample/marker overlap and maybe complex join.
			# need to check if overlapping marker sequences are identical
			cat("Complex join\n")
			if(!identical(master@markers[markerOverlap],resultList[[i]]@markers[markerOverlap]))  {
				#STOP
				stop("Cannot combine: markers with same name have different sequences.\n")
			}
			master <- mergeMlgtResults.complex(master, resultList[[i]])
		} else {					# no overlap in markers, 'simple' join
			cat("Simple join\n")
			master <- mergeMlgtResults.simple(master, resultList[[i]])
		}	
	}
	return(master)
}


#' @name my.mlgt.Result
#' @title An example \code{\link{mlgtResult}} object.
#' @description This is the result of running \code{\link{mlgt}} on the sample data given in the README.
#' @docType data
#' @usage my.mlgt.Result
#' @format a \code{\link{mlgtResult}} object.
#' @source Package mlgt
#' @author Dave T. Gerrard, 2012-04-01
NULL






