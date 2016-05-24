#' Methods for Editing or Converting Output from simFossilRecord
#'
#' These are a set of functions available for manipulating, translating
#' and editing the objects of class \code{fossilRecordSimulation} output
#' from function \code{simFossilRecord}.

#' @name simFossilRecordMethods

#' @details
#' These functions exist to manipulate \code{fossilRecordSimulation} objects
#' output from \code{simFossilRecord}, particularly so that they can be interfaced
#' with functions in library \code{paleotree} in the same way that output from the
#' deprecated 'legacy' simulation function \code{simFossilTaxa} was used.
#'
#' \code{timeSliceFossilRecord} takes a given \code{fossilRecordSimulation} object
#' and 'slices' the data to remove any events that occur after the given
#' \code{sliceTime} and make it so any taxa still alive as of \code{sliceTime}
#' are now listed as extant.
#'
#' \code{fossilRecord2fossilTaxa} converts a \code{fossilRecordSimulation} object
#' to the flat table format of taxon data as was originally output by deprecated function 
#' \code{simFossilTaxa}, and can be taken as input by a number of \code{paleotree} functions such as
#' \code{sampleRanges},\code{taxa2phylo} and \code{taxa2cladogram}.
#'
#' \code{fossilRecord2fossilRanges} converts a \code{fossilRecordSimulation} object
#' to the flat table format of observed taxon ranges, as is typically output by processing
#' \code{simFossilRecord} simulation output with \code{paleotree} function
#' \code{sampleRanges}.
#'

#' @param fossilRecord A list object output by \code{simFossilRecord}, often composed
#' of multiple elements, each of which is data for 'one taxon', with the first
#' element being a distinctive six-element vector composed of numbers, corresponding
#' to the six fields in tables output by the deprecated function \code{simFossilTaxa}.

#' @param sliceTime The date to slice the \code{simFossilRecord} output at, given
#' in time-units before the modern, on the same scale as the input \code{fossilRecord}.

#' @param shiftRoot4TimeSlice Should the dating of events be shifted, so that the
#' date given for \code{sliceTime} is now 0, or should the dates not be shifted,
#' so that they remain on the same scale as the input? This argument accepts a
#' logical TRUE or FALSE, but also accepts the string \code{"withExtantOnly"},
#' which will only 'shift' the time-scale if living taxa are present, as
#' determined by having ranges that overlap within \code{tolerance} of \code{sliceTime}.

#' @param tolerance A small number which sets a range around the \code{sliceTime} within
#' which taxa will be considered extant for the purposes of output.

#' @param modern.samp.prob The probability that a taxon is sampled at the modern time
#' (or, for \code{timeSliceFossilRecord}, the time at which the simulation data is
#' slice). Must be a number between 0 and 1. If 1, all taxa that survive to the modern
#' day (to the \code{sliceTime}) are sampled, if 0, none are.

#' @param merge.cryptic If \code{TRUE}, sampling events for cryptic taxon-units (i.e.
#' those in the same cryptic complex) will be merged into sampling events for a single
#' taxon-unit (with the name of the first taxon in that cryptic complex).
 
#' @param ranges.only If \code{TRUE} (the default), \code{fossilRecord2fossilRanges}
#' will return the dates of the first and last sampled occurrences of each taxon-unit
#' (i.e. the stratigraphic range of each taxon). If \code{FALSE}, instead a list will be
#' output, with each element being a vector of dates for all sampling events of each taxon-unit.

#' @return
#' Depends on the function and the arguments given. See Details.

#' @aliases timeSliceFossilRecord fossilRecord2fossilTaxa fossilRecord2fossilRanges 

#' @seealso
#' \code{\link{simFossilRecord}}

#' @author 
#' David W. Bapst.

#' @examples
#' 
#' set.seed(44)
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1,
#' 	nTotalTaxa=c(20,30) ,nExtant=0, plot=TRUE)
#' 
#' # time-slicing
#' 
#' # let's try slicing this record at 940 time-units
#' slicedRecord<-timeSliceFossilRecord(fossilRecord = record, sliceTime = 940)
#' # and let's plot it
#' divCurveFossilRecordSim(slicedRecord)
#' 
#' # now with shiftRoot4TimeSlice = TRUE to shift the root age
#' slicedRecord<-timeSliceFossilRecord(fossilRecord = record, sliceTime = 940,
#' 	shiftRoot4TimeSlice = TRUE)
#' # and let's plot it
#' divCurveFossilRecordSim(slicedRecord)
#' 
#' # plot look a little different due to how axis limits are treated...
#' # notice that in both, 'modern' (extant) taxa are sampled with probability = 1
#' 	#let's try it again, make that probability = 0
#' 
#' # now with shiftRoot4TimeSlice=TRUE
#' slicedRecord<-timeSliceFossilRecord(fossilRecord = record, sliceTime = 940,
#' 	shiftRoot4TimeSlice = TRUE, modern.samp.prob = 0)
#' # and let's plot it
#' divCurveFossilRecordSim(slicedRecord)
#' 
#' ############################
#' 
#' # converting to taxa objects and observed ranges
#' 
#' # convert to taxa data
#' taxa<-fossilRecord2fossilTaxa(record)
#' # convert to ranges
#' ranges<-fossilRecord2fossilRanges(record)
#' 
#' # plot diversity curves with multiDiv
#' multiDiv(list(taxa,ranges),plotMultCurves=TRUE)
#' # should look a lot like what we got earlier
#' 
#' # get the cladogram we'd obtain for these taxa with taxa2cladogram
#' cladogram<-taxa2cladogram(taxa,plot=TRUE)
#' 
#' # now get the time-scaled phylogenies with taxa2phylo
#' 
#' # first, with tips extending to the true times of extinction
#' treeExt<-taxa2phylo(taxa,plot=TRUE)
#' 
#' # now, with tips extending to the first appearance dates (FADs) of taxa
#' 	# get the FADs from the ranges
#' FADs<-ranges[,1]
#' treeFAD<-taxa2phylo(taxa,FADs,plot=TRUE)
#' 
#' @rdname simFossilRecordMethods
#' @export
timeSliceFossilRecord<-function(fossilRecord, sliceTime, shiftRoot4TimeSlice=FALSE,
		modern.samp.prob=1, tolerance=10^-4){
	#take a fossilRecord data object and cut it at some specific date
	#
	# CHECKS
	checkResult<-checkFossilRecord(fossilRecord)
	#check shiftRoot4TimeSlice
	shiftPar<-c(TRUE,FALSE,"withExtantOnly")
	shiftRoot4TimeSlice<-shiftPar[pmatch(shiftRoot4TimeSlice,shiftPar)]
	if(is.na(shiftRoot4TimeSlice)){
		stop("shiftRoot4TimeSlice must be a logical or the string 'withExtantOnly'")}
	#
	#drop all taxa that originate after the sliceTime
	droppers<-sapply(fossilRecord,function(x) x[[1]][3]<sliceTime)
	fossilRecord<-fossilRecord[!droppers]
	#
	#remove all sampling events after sliceTime
		for(i in 1:length(fossilRecord)){
			#remove all sampling events after sliceTime
			fossilRecord[[i]][[2]]<-fossilRecord[[i]][[2]][fossilRecord[[i]][[2]]>=sliceTime]
		}
	# adjusting time, making taxa extant
		# need to first test if there are extant taxa or not
	isAlive<-sapply(fossilRecord,function(x){
		if(is.na(x[[1]][4])){
			TRUE
		}else{
			(sliceTime-x[[1]][4])>tolerance
		}})
	#browser()
	if(shiftRoot4TimeSlice=="withExtantOnly"){
		if(any(isAlive)){
			shiftRoot4TimeSlice<-TRUE
		}else{
			shiftRoot4TimeSlice<-FALSE
			}
		}
	#
	# if shiftRoot4TimeSlice, then the whole thing shifts so time=0 is slice time
	if(shiftRoot4TimeSlice){	
		#adjust all dates so cutdate becomes 0
		for(i in 1:length(fossilRecord)){
			#adjust all dates so cutdate becomes 0
				#if stillAlive, replace 4:5 with 0,1
			if(isAlive[i]){
				#turn all taxa that went extinct after sliceTime so they are still alive
				fossilRecord[[i]][[1]][3]<-fossilRecord[[i]][[1]][3]-sliceTime
				fossilRecord[[i]][[1]][4:5]<-c(0,1)
			}else{
				fossilRecord[[i]][[1]][3:4]<-fossilRecord[[i]][[1]][3:4]-sliceTime
				}
			fossilRecord[[i]][[2]]<-fossilRecord[[i]][[2]]-sliceTime
			}
		modernTime<-0
	# if shiftRoot4TimeSlice=FALSE, then simply replace all extant taxa with
		# LADs at sliceTime and score as extant
	}else{
		for(i in 1:length(fossilRecord)){
			if(isAlive[i]){
				fossilRecord[[i]][[1]][4:5]<-c(sliceTime,1)
				}
			}
		modernTime<-sliceTime
		}
	#
	# sample at modern based on modern.samp.prob
	whichExtant<-which(sapply(fossilRecord,function(x) x[[1]][5]==1))
	nLive<-length(whichExtant)
	liveSampled<-as.logical(rbinom(n=nLive, size=1, prob=modern.samp.prob))
	whichSampled<-whichExtant[liveSampled]
	#
	#add sampling event at modern
	for(i in whichSampled){
		fossilRecord[[i]][[2]]<-c(fossilRecord[[i]][[2]],modernTime)
		}
	#
	# make sure it has the right class
	class(fossilRecord)<-'fossilRecordSimulation'
	# 
	return(fossilRecord)
	}

#' @rdname simFossilRecordMethods
#' @export
fossilRecord2fossilTaxa<-function(fossilRecord){
	# CHECKS
	checkResult<-checkFossilRecord(fossilRecord)
	# a function that transforms a simfossilrecord to a taxa object
	taxaConvert<-t(sapply(fossilRecord,function(x) x[[1]]))	
	rownames(taxaConvert)<-names(fossilRecord)
	return(taxaConvert)
	}

#' @rdname simFossilRecordMethods
#' @export	
fossilRecord2fossilRanges<-function(fossilRecord, merge.cryptic=TRUE, ranges.only = TRUE){
	# a function that transforms a simfossilrecord to a set of ranges (like from sampleRanges)
		# merge.cryptic = TRUE or FALSE
		# ranges.only or sampling times?
	# CHECKS
	# browser()
	checkResult<-checkFossilRecord(fossilRecord)
	#
	sampData<-lapply(fossilRecord,function(x) x[[2]]) 
	#get sampOcc : separate out the sampling events
	sampOcc<-lapply(fossilRecord,function(x) x[[2]])
	names(sampOcc)<-names(fossilRecord)
	#merge cryptic taxa
	if(merge.cryptic){
		taxonIDs<-sapply(fossilRecord,function(x) x[[1]][1])
		cryptIDs<-sapply(fossilRecord,function(x) x[[1]][6])
		for(i in 1:length(fossilRecord)){
			if(taxonIDs[i]==cryptIDs[i]){
				#if its the original taxon, collect all sampling events
					# for this cryptic complex into one pool
				# browser()
				sampOcc[[i]]<-unlist(c(sampOcc[taxonIDs[i]==cryptIDs]))
				#check that its a vector
				if(is.list(sampOcc[[i]])){
					stop("sampling data for taxa is not coercing correctly to a vector")}
			}else{
				#if its a cryptic taxon that didn't found the complex, erase its data
				sampOcc[[i]]<-NA
				}
			}
		}
	sampOcc[sapply(sampOcc,length)==0]<-NA
	#convert sampling events to FADs and LADs
	if(ranges.only){
		ranges<-cbind(sapply(sampOcc,max),sapply(sampOcc,min))
		rownames(ranges)<-names(sampOcc)
		colnames(ranges)<-c("FAD","LAD")
		result<-ranges
	}else{
		result<-sampOcc
		}
	return(result)
	}

	
# don't export
checkFossilRecord<-function(fossilRecord){
	if(!inherits(fossilRecord,"fossilRecordSimulation")){
		stop("fossilRecord object is not of class 'fossilRecordSimulation'")}
	if(any(sapply(fossilRecord,length)!=2)){
		stop("fossilRecord object has taxon entries with more or less than two elements")}
	if(any(sapply(fossilRecord,function(x) length(x[[1]]))!=6)){
		stop("fossilRecord object has taxon entries with more or less than six elements in first element")}		
	return(TRUE)
	}

