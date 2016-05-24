#' Construct a Stochastic Sequenced Time-List from an Unsequenced Time-List
#'
#' This function randomly samples from a timeList object (i.e. a list composed of a matrix of interval start and end
#' dates and a matrix of taxon first and last intervals), to find a set of taxa and intervals that do not overlap,
#' output as a new timeList object.

#' @details
#' Many analyses of diversification and sampling in the fossil record require a dataset composed of sequential non-overlappling intervals,
#' but the nature of the geologic record often makes this difficult, with taxa from different regions, environments and sedimentary basins
#' having first and last appearances placed in entirely in-congruent systems of chronostratigraphic intervals. While one option is to convert
#' such occurrences to a single, global stratigraphic system, this may still result in overlapping intervals when fossil collections are poorly
#' constrained stratigraphically. (For example, this may often be the case in global datasets.) 

#' This function offers an approach to avoid this issue in large datasets by randomly subsampling
#' the available taxa and intervals to produce stochastic
#' sets of ranges composed of data drawn from non-overlapping intervals. 
#'
#' This function is stochastic and thus should be set for many runs to produce many such solutions. Additionally,
#' all solutions found are returned, and users may wish to sort amongst these to maximize the number of intervals and 
#' number of taxa returned. A single solution which maximizes returned taxa and intervals may not be a precise enough approach
#' to estimating sampling rates, however, given the uncertainty in data. Thus, many runs should always be considered.
#'
#' By default, solutions are searched for without consideration to the length of intervals used (i.e. the selection of intervals is 'unweighted').
#' Alternatively, we can 'weight' selection toward the smallest intervals in the set, using the argument \code{weightSampling}. Smaller
#' intervals presumably overlap less and thus should retain more taxa and intervals of more equal length. However, in practise with empirical datasets,
#' the package author finds these approaches do not seem to produce very different estimates.
#'
#' For some datasets, many solutions found using seqTimeList may return infinite sampling values. This is often due to saving too many taxa
#' found in single intervals to the exclusion of longer-ranging taxa (see the example). This excess of single interval taxa is a clear artifact
#' of the randomized seqTimeList procedure and such solutions should probably be ignored.
 
#' @param timeList A list composed of two matrices, giving interval start and end 
#' dates and taxon first and last occurrences within those intervals. Some intervals
#' are expected to overlap (thus necessitating the use of this function), and datasets
#' lacking overlapping intervals will return an error message.

#' @param nruns Number of new timeList composed of non-overlapping intervals produced.

#' @param weightSampling If TRUE, weight sampling of new intervals toward smaller intervals. FALSE by default.

#' @return
#' A list, composed of three elements: \code{nIntervals} which is a vector of the
#' number of intervals in each solution, \code{nTaxa} which is a vector of the number of
#' taxa in each solution and \code{timeLists} which is a list composed of each new
#' timeList object as an element.

#' @seealso 
#' Resulting time-lists can be analyzed with \code{\link{freqRat}}, \code{\link{durationFreq}}, etc.
#' 
#' Additionally, \code{\link{binTimeData}} can be useful for simulating interval data.

#' @author David W. Bapst

#' @examples
#' #Simulate some fossil ranges with simFossilRecord
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(60,80), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.1)
#' 
#' #Now let's use binTimeData to get ranges in discrete overlapping intervals
#'     #via pre-set intervals input
#' presetIntervals <- cbind(c(1000,995,990,980,970,975,960,950,940,930,900,890,888,879,875),
#'   c(995,989,960,975,960,950,930,930,930,900,895,888,880,875,870))
#' rangesDisc1 <- binTimeData(rangesCont,int.times=presetIntervals)
#' 
#' seqLists<-seqTimeList(rangesDisc1,nruns=10)
#' seqLists$nTaxa
#' seqLists$nIntervals
#' 
#' #apply freqRat as an example analysis
#' sapply(seqLists$timeLists,freqRat)
#' 
#' #notice the zero and infinite freqRat estimates? What's going on?
#' 
#' freqRat(seqLists$timeLists[[4]],plot=TRUE)
#' 
#' #too few taxa of two or three interval durations for the ratio to work properly
#'     #perhaps ignore these estimates
#' 
#' #with weighted selection of intervals
#' seqLists<-seqTimeList(rangesDisc1,nruns=10,weightSampling=TRUE)
#' 
#' seqLists$nTaxa
#' seqLists$nIntervals
#' sapply(seqLists$timeLists,freqRat)
#' 
#' #didn't have much effect in this simulated example

#' @export seqTimeList
seqTimeList<-function(timeList,nruns=100,weightSampling=FALSE){
	timeList[[1]]<-as.matrix(timeList[[1]])
	timeList[[2]]<-as.matrix(timeList[[2]])
	#let's try to get sampling rates out of a timeList with overlapping intervals
	#test for overlap
	lap<-any(!apply(timeList[[1]],1,function(x)
		2>sum(x[1]>(timeList[[1]][,2]) & x[2]<(timeList[[1]][,1]))))
	if(!lap){stop("Your timeList doesn't have any overlap??")}
	#let's make a list of valid intervals
		#first, all non-overlapping intervals
	noLap<-which(apply(timeList[[1]],1,function(x) 
		2>sum(x[1]>(timeList[[1]][,2]) & x[2]<(timeList[[1]][,1]))))
	#now let's search iteratively for interval solutions
	nInts<-nTaxa<-numeric()
	savedLists<-list(1:nruns)
	for(i in 1:nruns){
		#set it up!
		valid<-noLap
		invalid<-numeric()
		if(length(valid)>0){
			candidate<-(1:nrow(timeList[[1]]))[-valid]
			}else{
			candidate<-(1:nrow(timeList[[1]]))
			}
		#LOOP
		while(length(candidate)>1){
			#next, let's pull out a random candidate interval
			if(weightSampling){
				weights<-order(apply(timeList[[1]][candidate,,drop=FALSE],1,diff))
				weights<-weights/sum(weights)
				draw<-sample(1:length(candidate),1,prob=weights)
			}else{
				draw<-sample(1:length(candidate),1)
				}
			#remove from candidate, add to valid
			valid<-c(valid,candidate[draw])		
			drawInt<-timeList[[1]][candidate[draw],]
			candidate<-candidate[-draw]
			#add all overlapping to invalid and remove from candidate
			candInt<-timeList[[1]][candidate,,drop=FALSE]
			badInt<-apply(candInt,1,function(x) drawInt[1]>x[2] & drawInt[2]<x[1])	
			invalid<-c(invalid,candidate[badInt])
			candidate<-candidate[!badInt]
			#
			#newInts<-timeList[[1]][valid,]
			#if(any(!apply(newInts,1,function(x)
			#	2>sum(x[1]>(newInts[,2]) & x[2]<(newInts[,1]))))){
			#		stop("WHAT!")}
			}
		if(length(candidate)==1){valid<-c(valid,candidate)}
		valid<-sort(valid)
		newInts<-timeList[[1]][valid,]
		# test that there be no overlaps
		if(any(!apply(newInts,1,function(x)	
			2>sum(x[1]>(newInts[,2]) & x[2]<(newInts[,1]))))){stop("WHAT!")}
		#
		validTaxa<-apply(timeList[[2]],1,function(x) any(x[1]==valid) & any(x[2]==valid))
		cropTaxa<-timeList[[2]][validTaxa,]
		newTaxa<-cbind(sapply(cropTaxa[,1],function(x) which(x==valid)),
		sapply(cropTaxa[,2],function(x) which(x==valid)))
		nInts[i]<-nrow(newInts)
		nTaxa[i]<-nrow(newTaxa)
		savedLists[[i]]<-list(intTimes=newInts,taxonTimes=newTaxa)
		}
	results<-list(nIntervals=nInts,nTaxa=nTaxa,timeLists=savedLists)
	return(results)
	}