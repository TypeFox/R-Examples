#' perCapitaRates
#' 
#' Calculates and plots per-capita origination and extinction rates from
#' sequential discrete-time taxon ranges, following Foote (2000).
#' 
#' @details
#' This function calculates the per-capita rates of taxonomic origination
#' and extinction from paleontological range data, as described by Foote 
#' (2000). These values are the instantaneous rate of either type of event
#' occurring per lineage time-units. Although Foote (2001) also presents a
#' number of alternative rates collected from the prior literature such as
#' the VanValen rate metrics, these are not implemented here, but could be
#' estimated using the matrix invisibly output by this function.
#'
#' The timeList object should be a list composed of two matrices, the first
#' matrix giving by-interval start and end times (in absolute time), the second
#' matrix giving the by-taxon first and last appearances in the intervals
#' defined in the first matrix, numbered as the rows. Absolute time should be
#' decreasing, while the intervals should be numbered so that the number
#' increases with time. Taxa alive in the modern should be either (a) listed 
#' in isExtant or (b) listed as last occurring in a time interval that 
#' begins at time 0 and ends at time 0. See the documentation for the time-scaling 
#' function \code{\link{bin_timePaleoPhy}} and the simulation function 
#' \code{\link{binTimeData}} for more information on formatting.
#'
#' Unlike some functions in paleotree, such as the diversity curve functions,
#' intervals must be both sequential and non-overlapping. The diversity curve
#' functions deal with such issues by assuming taxa occur from the base of the
#' interval they are first found in until the end of the last interval they
#' are occur in. This inflation of boundary crossers could badly bias estimates
#' of per-capita diversification rates.

#' @inheritParams DiversityCurves

#' @param plot If true, the per-capita origination and extinctions rates are 
#' plotted for each interval. Rates which cannot be calculated for an interval
#' will not be plotted, thus appearing as a gap in the plotted graph. The
#' author takes no responsibility for the aesthetics of this plot.

#' @param logRates If true, rates plotted on log scale.

#' @param drop.extant Drops all extant taxa from a dataset before
#' calculating per-capita origination and extinction rates.

#' @param isExtant A vector of TRUE and FALSE values, same length as the
#' number of taxa in the second matrix of timeList, where TRUE values indicate
#' taxa that are alive in the modern day (and thus are boundary crossers which
#' leave the most recent interval). By default, this argument is NULL and instead
#' which taxa are extant is inferred based on which taxa occur in an interval
#' with start and end times both equal to zero. See details.

#' @param jitter If TRUE (default) the extinction rate will be plotted slightly
#' ahead of the origination rate on the time axis, so the two can be differentiated.

#' @param legendPosition The position of a legend indicating which line is
#' origination rate and which is extinction rate on the resulting plot. This
#' is given as the possible positions for argument 'x' of the function 
#'\code{\link{legend}}, and by default is "topleft", which will be generally
#' useful if origination and extinction rates are initially low. If 
#' legendPosition is NA, then a legend will not be plotted.

#' @return This function will invisibly return a ten column matrix,
#' where the number of rows is equal to the number of intervals. The
#' first two columns are interval start and end times and the third
#' column is interval length. The fourth through eighth column is the
#' four fundamental classes of taxa from Foote (2001): Nbt, NbL, NFt,
#' NFL and their sum, N. The final two columns are the per-capita
#' rates estimated for each interval in units per lineage time-units;
#' the ninth column is the origination rate ('pRate') and the tenth
#' column is the extinction rate ('qRate').

#' @seealso \code{\link{DiversityCurves}}, \code{\link{SamplingConv}}

#' @references
#' Foote, M. 2000 Origination and extinction components of taxonomic diversity:
#' general problems. Pp. 74--102. In D. H. Erwin, and S. L. Wing, eds. Deep
#' Time: Paleobiology's Perspective. The Paleontological Society, Lawrence,
#' Kansas.

#' @examples
#' 
#' #with the retiolinae dataset
#' data(retiolitinae)
#' perCapitaRates(retioRanges)
#'
#' #Simulate some fossil ranges with simFossilRecord
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(80,100), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #Now let's use binTimeData() to bin in intervals of 5 time units
#' rangesDisc <- binTimeData(rangesCont,int.length=5)
#' #and get the per-capita rates
#' perCapitaRates(rangesDisc)
#' #on a log scale
#' perCapitaRates(rangesDisc,logRates=TRUE)
#' 
#' #get mean and median per-capita rates
#' res<-perCapitaRates(rangesDisc,plot=FALSE)
#' apply(res[,c("pRate","qRate")],2,mean,na.rm=TRUE)
#' apply(res[,c("pRate","qRate")],2,median,na.rm=TRUE)
#' 
#' #with modern taxa
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(80,100))
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.5,,modern.samp.prob=1)
#' #Now let's use binTimeData() to bin in intervals of 5 time units
#' rangesDisc <- binTimeData(rangesCont,int.length=5)
#' #and now get per-capita rates
#' perCapitaRates(rangesDisc)
#'
#' @export
perCapitaRates<-function(timeList,plot=TRUE,logRates=FALSE,drop.extant=FALSE,isExtant=NULL,jitter=TRUE,legendPosition="topleft"){
	#
	#this function estimates per-capita rates for binned intervals from discrete interval range data
		#based on Foote, 2000
	#input is a list with (1) interval times matrix and (2) species FOs and LOs
	#time interval starts and ends can be pre-input as a 2 column matrix
		#HOWEVER this could be pretty misleading!
		#standing richness may never be high as the apparent richness of some bins
	#if there is a single interval that is start=0 and end=0, it is dropped and all appearances in this interval are either
		#(a) melded onto the next most recent interval
		#(b) dropped if drop.extant=TRUE
	#output is a matrix of int-start, int-end, p, q
	#example
		#timeList<-rangesDisc;int.times=NULL;plot=TRUE;plotLogRates=FALSE;timelims=NULL
			#drop.extant=FALSE;isExtant=NULL;split.int=TRUE
	#
	if(!inherits(timeList[[1]],"matrix")){
		if(inherits(timeList[[1]],"data.frame")){
			timeList[[1]]<-as.matrix(timeList[[1]])
		}else{
			stop("timeList[[1]] not of matrix or data.frame format")
			}
		}
	if(!inherits(timeList[[2]],"matrix")){
		if(inherits(timeList[[2]],"data.frame")){
			timeList[[2]]<-as.matrix(timeList[[2]])
		}else{
			stop("timeList[[2]] not of matrix or data.frame format")
			}
		}
	intMat<-timeList[[1]]	#the intervals the DATA is given in
	timeData<-timeList[[2]]
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(intMat,1,diff)>0)){stop("timeList[[1]] not in intervals in time relative to modern")}
	if(any(intMat[,2]<0)){stop("Some dates in timeList[[1]] <0 ?")}
	if(any(apply(timeData,1,diff)<0)){stop("timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Some dates in timeList[[2]] <0 ?")}
	#get rid of modern intervals
	modInt<-which(intMat[,1]==0)
	if(length(modInt)>0){
		if(drop.extant){
			timeData<-timeData[!(timeData[,2]==modInt),]
		}else{
			if(is.null(isExtant)){isExtant<-timeData[,2]==modInt}
			altMod<-which(!intMat[,1]==0 & intMat[,2]==0)
			timeData[timeData[,1]==modInt,1]<-altMod
			timeData[timeData[,2]==modInt,2]<-altMod
			}
		intMat<-intMat[-modInt,]
		}
	if(is.null(isExtant)){isExtant<-rep(FALSE,nrow(timeData))}	
	# make sure data is absolutely sequential and overlapping
	isSequential_1<-all(apply(intMat,2,function(x) identical(x,rev(sort(x)))))
	isNonRep<-all(apply(intMat,2,function(x) identical(length(x),length(unique(x)))))
	isSequential_2<-all(sapply(2:nrow(intMat),function(x) intMat[x,1]==intMat[x-1,2]))
	if(!(isSequential_1 & isSequential_2 & isNonRep)){
		stop("Sorry, intervals need to be sequential, ordered and non-overlapping.")}
	#now get int length to modify rates with
	intlen<-(-apply(intMat,1,diff))
	#modify timeData if there are any extant taxa
	if(any(isExtant)){timeData[isExtant,2]<-nrow(intMat)+1}
	#get the four values
	Nbt<-sapply(1:nrow(intMat),function(x) sum(timeData[,1]<x & timeData[,2]>x))
	NbL<-sapply(1:nrow(intMat),function(x) sum(timeData[,1]<x & timeData[,2]==x))
	NFt<-sapply(1:nrow(intMat),function(x) sum(timeData[,1]==x & timeData[,2]>x))	
	NFL<-sapply(1:nrow(intMat),function(x) sum(timeData[,1]==x & timeData[,2]==x))
	N<-Nbt+NbL+NFt+NFL	#total diversity
	#now calculate rates
	pRate<-ifelse(Nbt==0,NA,-log(Nbt/(Nbt+NFt))/intlen)
	qRate<-ifelse(Nbt==0,NA,-log(Nbt/(Nbt+NbL))/intlen)
	if(plot){
		int.start<-intMat[,1];int.end<-intMat[,2]
		times1<-c(int.start,(int.end+((int.start-int.end)/100)))
		#add a jigger so rates don't overlap
		jigger<-min(diff(times1))/10
		p1<-c(pRate,pRate)[order(times1)]
		q1<-c(qRate,qRate)[order(times1)]
		times1<-sort(times1)
		if(logRates){
			p1[!(p1>0)]<-NA
			q1[!(q1>0)]<-NA
			ylims<-c(min(c(p1,q1),na.rm=TRUE),(max(c(p1,q1),na.rm=TRUE))*1.5)
			plot(times1,p1,type="l",log="y",col=4,lwd=2,lty=5,
				xlim=c(max(times1),max(0,min(times1))),ylim=ylims,
				xlab="Time (Before Present)",ylab="Instantaneous Per-Capita Rate (per Ltu)")
		}else{
			ylims<-c(min(c(p1,q1),na.rm=TRUE),(max(c(p1,q1),na.rm=TRUE))*1.2)
			plot(times1,p1,type="l",col=4,lwd=2,lty=5,
				xlim=c(max(times1),max(0,min(times1))),ylim=ylims,
				xlab="Time (Before Present)",ylab="Instantaneous Per-Capita Rate (per Ltu)")
			}
		if(jitter){
			lines(times1+jigger,q1,col=2,lwd=2,lty=2)
		}else{
			lines(times1,q1,col=2,lwd=2,lty=2)
			}
		if(!is.na(legendPosition)){legend(x=legendPosition,legend=c("Origination","Extinction"),lty=c(5,2),lwd=2,col=c(4,2))}
		}
	res<-cbind(intMat,intlen,Nbt,NbL,NFt,NFL,N,pRate,qRate)
	return(invisible(res))
	}