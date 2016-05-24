#' Bin Simulated Temporal Ranges in Discrete Intervals
#' 
#' Converts a matrix of simulated continuous-time first occurrences and last
#' occurrences for fossil taxa into first and last occurrences given in some set
#' of simulated/placed discrete-time intervals, which is output along with
#' information of the dates of the given intervals.
#' 
#' @details This function takes a simulated matrix of per-taxon first and last
#' occurrences and, by dividing the time-scale into time intervals of non-zero
#' length, lists taxon occurrences within those interval. By default, a set of
#' sequential non-overlapping time-interval of equal non-zero length are used,
#' with the length controlled by the argument int.length.
#' 
#' Alternatively, a two column matrix of interval start and end times to be
#' used can be input via the argument int.times. None of these intervals can
#' have a duration (temporal length) greater than zero. If a first or last
#' appearance in the input range data could fit into multiple intervals (i.e.
#' the input discrete time intervals are overlapping), then the appearance data
#' is placed in the interval of the shortest duration. When output, the
#' interval times matrix (see below) will be sorted from first to last.
#' 
#' This function is SPECIFICALLY for simulating the effect of having a discrete
#' time-scale for analyses using simulations. This function should not be used
#' for non-simulations uses, such as binning temporal occurrences for analyses
#' of real data. In those case, the temporal ranges (which, in real data, will
#' probably be given as discrete time intervals) should already be tabulated
#' within discrete intervals prior to use in paleotree. The user should place
#' the temporal information in a list object, as described for the output of
#' binTimeData (see below).
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero. However, the numbering of
#' intervals giving in the output increases with time, as these are numbered
#' relative to each other, from first to last.
#' 
#' As of version 1.7 of paleotree, taxa which are extant as indicated in timeData as being
#' in a time interval bounded (0,0), unless time-bins are preset using
#' argument int.times (prior to version 1.5 they were erroneously listed as
#' NA).
#' 
#' @param timeData Two-column matrix of simulated first and last occurrences in
#' absolute continuous time

#' @param int.length Time interval length, default is 1 time unit

#' @param start Starting time for calculating the intervals.

#' @param int.times A two column matrix with the start and end times of the
#' intervals to be used.

#' @return A list containing: \item{int.times}{A 2 column matrix with the start
#' and end times of the intervals used; time decreases relative to the
#' present.} \item{taxon.times}{A 2 column matrix with the first and last
#' occurrences of taxa in the intervals listed in int.times, with numbers
#' referring to the row of int.times.}

#' @seealso \code{\link{simFossilRecord}}, \code{\link{sampleRanges}},
#' \code{\link{taxicDivCont}}

#' @author David W. Bapst

#' @examples
#' 
#' #Simulate some fossil ranges with simFossilRecord
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #Now let's use binTimeData() to bin in intervals of 1 time unit
#' rangesDisc <- binTimeData(rangesCont,int.length=1)
#' #plot with taxicDivDisc()
#' equalDiscInt <- taxicDivDisc(rangesDisc)
#' 
#' #example with pre-set intervals input (including overlapping)
#' presetIntervals <- cbind(c(1000,990,970,940),c(980,970,950,930))
#' rangesDisc1 <- binTimeData(rangesCont,int.times=presetIntervals)
#' taxicDivDisc(rangesDisc1)
#' #now let's plot diversity with (different) equal length intervals used above
#' taxicDivDisc(rangesDisc1,int.times=equalDiscInt[,1:2])
#' 
#' #example with extant taxa
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40))
#' taxa<-fossilRecord2fossilTaxa(record)
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.5,,modern.samp.prob=1)
#' #Now let's use binTimeData() to bin in intervals of 1 time unit
#' rangesDisc <- binTimeData(rangesCont,int.length=1)
#' #plot with taxicDivDisc()
#' taxicDivDisc(rangesDisc)
#' 
#' #example with pre-set intervals input (including overlapping)
#' presetIntervals <- cbind(c(40,30,20,10),c(30,20,10,0))
#' rangesDisc1 <- binTimeData(rangesCont,int.times=presetIntervals)
#' taxicDivDisc(rangesDisc1)
#' 
#' @export binTimeData
binTimeData<-function(timeData,int.length=1,start=NA,int.times=NULL){
	#bin temporal data
	#input: continuous time data (two column of FADs and LADs)
	#output: a list with two 2-col matrices as elements, bin-times and taxon occurences
			#intervals, UNLIKE TIME, always go up (earliest is 1 and increase...)
		#arbitrarily starts bin at the first fad; this can be changed by setting 'start'
			#start must be greater than max(timeData)
			#the last bin is cut off at zero (present day)
	#x<-c(0,runif(99));timeData<-cbind(x+rexp(100),x);int.length=1;start=NA;int.times=NULL
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data?")}
	if(any(timeData[,1]<timeData[,2])){stop("timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Some dates in timeData <0 ?")}
	if(is.null(int.times)){
		if(is.na(start)){start<-max(timeData)+int.length}else{if(start<max(timeData)){stop("Error:Start<max(timeData)?")}}
		end<-start-(ceiling((start-min(timeData))/int.length)+1)*int.length
		bins<-seq(start,end,by=-int.length)
		bins<-unique(ifelse(bins<0,0,bins))	#get rid of any extra zeroes or negative numbers
		fads<-sapply(timeData[,1],function(x) which(bins<x)[1]-1)
		lads<-sapply(timeData[,2],function(x) which(bins<x)[1]-1)
		if(any(timeData[,1]==0) | any(timeData[,2]==0)){
			bins<-c(bins,0)
			fads[timeData[,1]==0]<-length(bins)-1
			lads[timeData[,2]==0]<-length(bins)-1
			}
		res<-list(int.times=cbind(int.start=bins[1:(length(bins)-1)],int.end=bins[2:length(bins)]),
			taxon.times=cbind(first.int=fads,last.int=lads))
	}else{
		int.durs<-int.times[,1]-int.times[,2]
		if(any(int.durs<=0)){stop("Some input time intervals have zero or negative durations?")}
		int.times<-int.times[order(int.durs),]
		Fint<-sapply(timeData[,1],function(x) which(apply(int.times,1,function(y) y[1]>=x & y[2]<x))[1])
		Lint<-sapply(timeData[,2],function(x) which(apply(int.times,1,function(y) y[1]>=x & y[2]<x))[1])
		if(any(int.times[,2]==0)){
			Fint[timeData[,1]==0]<-which(int.times[,2]==0)[1]
			Lint[timeData[,2]==0]<-which(int.times[,2]==0)[1]
			}
		taxon.times<-cbind(first.int=Fint,last.int=Lint)
		rownames(taxon.times)<-rownames(timeData)
		taxon.times<-taxon.times[!apply(taxon.times,1,function(x) any(is.na(x))),]
		new.order<-rank(-int.times[,1])	
		taxon.times[,1]<-new.order[taxon.times[,1]]
		taxon.times[,2]<-new.order[taxon.times[,2]]
		int.times<-int.times[order(-int.times[,1]),]
		res<-list(int.times=int.times,taxon.times=taxon.times)
		}
	return(res)
	}
