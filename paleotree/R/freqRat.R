#' Frequency Ratio Method for Estimating Sampling Probability
#' 
#' Estimate per-interval sampling probability in the fossil record from a set
#' of discrete-interval taxon ranges using the frequency-ratio method described
#' by Foote and Raup (1996). Can also calculate extinction rate per interval from
#' the same data distribution.
#' 
#' @details This function uses the frequency ratio ("freqRat") method of Foote and Raup
#' (1996) to estimate the per-interval sampling rate for a set of taxa. This
#' method assumes that intervals are of fairly similar length and that
#' taxonomic extinction and sampling works similar to homogenous Poisson
#' processes. These analyses are ideally applied to data from single
#' stratigraphic section but potentially are applicable to regional or global
#' datasets, although the behavior of those datasets is less well understood.
#' 
#' The frequency ratio is a simple relationship between the number of taxa
#' observed only in a single time interval (also known as singletons), the
#' number of taxa observed only in two time intervals and the number of taxa
#' observed in three time intervals. These respective frequencies, respectively
#' f1, f2 and f3 can then be related to the per-interval sampling probability
#' with the following expression.
#' 
#' \eqn{Sampling.Probability = (f2^2)/(f1*f3)}
#' 
#' This frequency ratio is generally referred to as the 'freqRat' in
#' paleobiological literature.
#' 
#' It is relatively easy to visually test if range data fits expectation that
#' true taxon durations are exponentially distributed by plotting the
#' frequencies of the observed ranges on a log scale: data beyond the
#' 'singletons' category should have a linear slope, implying that durations
#' were originally exponentially distributed. (Note, a linear scale is used for
#' the plotting diagram made by this function when 'plot' is TRUE, so that plot
#' cannot be used for this purpose.)
#' 
#' The accuracy of this method tends to be poor at small interval length and
#' even relatively large sample sizes. A portion at the bottom of the examples
#' in the help file examine this issue in greater detail with simulations. This
#' package author recommends using the ML method developed in Foote (1997)
#' instead, which is usable via the function \code{\link{make_durationFreqDisc}}.
#' 
#' As extant taxa should not be included in a freqRat calculation, any taxa
#' listed as being in a bin with start time 0 and end time 0 (and thus being
#' extant without question) are dropped before the model fitting it performed.
#' 
#' @param timeData A 2 column matrix with the first and last occurrences of taxa
#' given in relative time intervals. If a list of length two is given for
#' timeData, such as would be expected if the output of binTimeData was
#' directly input, the second element is used.

#' @param calcExtinction If TRUE, the per-interval, per-lineage extinction rate 
#' is estimated as the negative slope of the log frequencies, ignoring single
#' hits (as described in Foote and Raup, 1996.)

#' @param plot If true, the histogram of observed taxon ranges is plotted, with
#' frequencies on a linear scale

#' @return This function returns the per-interval sampling probability as the
#' "freqRat", and estimates 

#' @author David W. Bapst

#' @seealso Model fitting methods in \code{\link{make_durationFreqDisc}} and \code{\link{make_durationFreqCont}}. 
#' Also see conversion methods in \code{\link{sProb2sRate}}, \code{\link{qsProb2Comp}}

#' @references Foote, M. 1997 Estimating Taxonomic Durations and Preservation
#' Probability. \emph{Paleobiology} \bold{23}(3):278--300.
#' 
#' Foote, M., and D. M. Raup. 1996 Fossil preservation and the stratigraphic
#' ranges of taxa. \emph{Paleobiology} \bold{22}(2):121--140.
#' @examples
#' 
#' #Simulate some fossil ranges with simFossilRecord
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #simulate a fossil record with imperfect sampling with sampleRanges
#' rangesCont <- sampleRanges(taxa,r=0.1)
#' #Now let's use binTimeData to bin in intervals of 5 time units
#' rangesDisc <- binTimeData(rangesCont,int.length=5)
#' 
#' #now, get an estimate of the sampling rate (we set it to 0.5 above)
#' #for discrete data we can estimate the sampling probability per interval (R)
#'     #i.e. this is not the same thing as the instantaneous sampling rate (r)
#' #can use sRate2sProb to see what we would expect
#' sRate2sProb(r=0.1,int.length=5)
#' #expect R = ~0.39
#' 
#' #now we can apply freqRat to get sampling probability
#' SampProb <- freqRat(rangesDisc,plot=TRUE)
#' SampProb
#'
#' #est. R = ~0.25 
#' #Not wildly accurate, is it?
#'
#' #can also calculate extinction rate per interval of time
#' freqRat(rangesDisc,calcExtinction=TRUE)
#'
#' #est. ext rate = ~0.44 per interval
#' #5 time-unit intervals, so ~0.44 / 5 = ~0.08 per time-unite
#' #That's pretty close to the generating value of 0.01, used in sampleRanges
#' 
#' \dontrun{
#' #################
#' #The following example code (which is not run by default) examines how 
#' 	#the freqRat estimates vary with sample size, interval length
#' 	#and compare it to using make_durationFreqDisc
#' 
#' #how good is the freqRat at 20 sampled taxa on avg?
#' set.seed(444)
#' r<-runif(100)
#' int.length=1
#' R<-sapply(r,sRate2sProb,int.length=1)	#estimate R from r, assuming stuff like p=q
#' ntaxa<-freqRats<-numeric()
#' for(i in 1:length(r)){
#' 	#assuming budding model
#' 	record<-simFossilRecord(p=0.1, q=0.1, r=r[i], nruns=1,
#' 		nSamp=c(15,25), nExtant=0, plot=TRUE)
#' 	ranges<-fossilRecord2fossilRanges(record)
#' 	timeList<-binTimeData(ranges,int.length=int.length)
#' 	ntaxa[i]<-nrow(timeList[[2]])
#' 	freqRats[i]<-freqRat(timeList)
#' 	}
#' plot(R,freqRats);abline(0,1)
#' #without the gigantic artifacts bigger than 1...
#' plot(R,freqRats,ylim=c(0,1));abline(0,1)
#' #very worrisome lookin'!
#' 
#' #how good is it at 100 sampled taxa on average?
#' set.seed(444)
#' r<-runif(100)
#' int.length=1
#' R<-sapply(r,sRate2sProb,int.length=1)
#' ntaxa<-freqRats<-numeric()
#' for(i in 1:length(r)){
#' 	#assuming budding model
#' 	record<-simFossilRecord(p=0.1, q=0.1, r=r[i], nruns=1,
#' 		nSamp=c(80,150), nExtant=0, plot=TRUE)
#' 	ranges<-fossilRecord2fossilRanges(record)
#' 	timeList<-binTimeData(ranges,int.length=int.length)
#' 	ntaxa[i]<-nrow(timeList[[2]])
#' 	freqRats[i]<-freqRat(timeList)
#' 	}
#' plot(R,freqRats,ylim=c(0,1));abline(0,1)
#' #not so hot, eh?
#' 
#' ################
#' #LETS CHANGE THE TIME BIN LENGTH!
#' 
#' #how good is it at 100 sampled taxa on average, with longer time bins?
#' set.seed(444)
#' r<-runif(100)
#' int.length<-10
#' R<-sapply(r,sRate2sProb,int.length=int.length)
#' ntaxa<-freqRats<-numeric()
#' for(i in 1:length(r)){
#' 	#assuming budding model
#' 	record<-simFossilRecord(p=0.1, q=0.1, r=r[i], nruns=1,
#' 		nSamp=c(80,150), nExtant=0, plot=TRUE)
#' 	ranges<-fossilRecord2fossilRanges(record)
#' 	timeList<-binTimeData(ranges,int.length=int.length)
#' 	ntaxa[i]<-nrow(timeList[[2]])
#' 	freqRats[i]<-freqRat(timeList)
#' 	}
#' plot(R,freqRats,ylim=c(0,1));abline(0,1)
#' #things get more accurate as interval length increases... odd, eh?
#' 
#' #how good is it at 20 sampled taxa on average, with longer time bins?
#' set.seed(444)
#' r<-runif(100)
#' int.length<-10
#' R<-sapply(r,sRate2sProb,int.length=int.length)
#' ntaxa<-freqRats<-numeric()
#' for(i in 1:length(r)){
#' 	#assuming budding model
#' 	record<-simFossilRecord(p=0.1, q=0.1, r=r[i], nruns=1,
#' 		nSamp=c(15,25), nExtant=0, plot=TRUE)
#' 	ranges<-fossilRecord2fossilRanges(record)
#' 	timeList<-binTimeData(ranges,int.length=int.length)
#' 	ntaxa[i]<-nrow(timeList[[2]])
#' 	freqRats[i]<-freqRat(timeList)
#' 	}
#' plot(R,freqRats,ylim=c(0,1));abline(0,1)
#' #still not so hot at low sample sizes, even with longer bins
#' 
#' ########################
#' #ML METHOD
#' 
#' #how good is the ML method at 20 taxa, 1 time-unit bins?
#' set.seed(444)
#' r<-runif(100)
#' int.length<-1
#' R<-sapply(r,sRate2sProb,int.length=int.length)
#' ntaxa<-ML_sampProb<-numeric()
#' for(i in 1:length(r)){
#' 	#assuming budding model
#' 	record<-simFossilRecord(p=0.1, q=0.1, r=r[i], nruns=1,
#' 		nSamp=c(15,25), nExtant=0, plot=TRUE)
#' 	ranges<-fossilRecord2fossilRanges(record)
#' 	timeList<-binTimeData(ranges,int.length=int.length)
#' 	ntaxa[i]<-nrow(timeList[[2]])
#'  likFun<-make_durationFreqDisc(timeList)
#'  ML_sampProb[i]<-optim(parInit(likFun),likFun,
#' 		lower=parLower(likFun),upper=parUpper(likFun),
#'      method="L-BFGS-B",control=list(maxit=1000000))[[1]][2]
#' 	}
#' plot(R,ML_sampProb);abline(0,1)
#' # Not so great due to likelihood surface ridges
#'  # but it returns values between 0-1
#' 
#' #how good is the ML method at 100 taxa, 1 time-unit bins?
#' set.seed(444)
#' r<-runif(100)
#' int.length<-1
#' R<-sapply(r,sRate2sProb,int.length=int.length)
#' ntaxa<-ML_sampProb<-numeric()
#' for(i in 1:length(r)){
#' 	#assuming budding model
#' 	record<-simFossilRecord(p=0.1, q=0.1, r=r[i], nruns=1,
#' 		nSamp=c(80,150), nExtant=0, plot=TRUE)
#' 	ranges<-fossilRecord2fossilRanges(record)
#' 	timeList<-binTimeData(ranges,int.length=int.length)
#' 	ntaxa[i]<-nrow(timeList[[2]])
#'  likFun<-make_durationFreqDisc(timeList)
#'  ML_sampProb[i]<-optim(parInit(likFun),likFun,
#' 		lower=parLower(likFun),upper=parUpper(likFun),
#'      method="L-BFGS-B",control=list(maxit=1000000))[[1]][2]
#' 	}
#' plot(R,ML_sampProb);abline(0,1)
#' #Oh, fairly nice, although still a biased uptick as R gets larger
#' 
#' }
#' 
#' @export freqRat
freqRat<-function(timeData,calcExtinction=FALSE,plot=FALSE){
	#timeData is discrete bin data, like from binTimeData
	if(length(timeData)==2){	#if a timeList matrix...
		timeData[[2]][(timeData[[1]][timeData[[2]][,2],2]==0),1]<-NA
		timeData<-timeData[[2]]
		}	
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(timeData,1,diff)<0)){
		stop("timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Some dates in timeList[[2]] <0 ?")}
	durations<-apply(timeData,1,diff)+1
	sumDur<-sapply(1:max(c(durations,3)),function(x) sum(durations==x))/length(durations)
	#get freqRat
	f1<-sumDur[1]
	f2<-sumDur[2]
	f3<-sumDur[3]
	freqRat<-(f2^2)/(f1*f3)
	names(freqRat)<-"freqRat"
	if(is.nan(freqRat)){
		message("Warning: Frequency distribution of input range data appears to violate model assumptions, producing a freqRat of zero over zero (NA)")
	}else{
		if(freqRat>1){
			message("Warning: Frequency distribution of input range data appears to violate model assumptions, producing an impossible freqRat greater than 1")
			}
		}
	#calculate extinction rate (rate of lineages going extinct per lineage, per interval)
	if(calcExtinction){
		logSD<-log(sumDur[-1])
		logSD[is.infinite(logSD)]<-NA
		intseq<-2:max(durations)
		reg<-lm(logSD~intseq)
		extRate<--coefficients(reg)[2]
		names(extRate)<-"extRate"
		freqRat<-c(freqRat,extRate)
		}
	#plotting
	if(plot){
		hist(durations,breaks=0:max(durations),xlab="Duration (time-units)",main="")
		}
	return(freqRat)
	}
