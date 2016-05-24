#' Sampling Taxon Ranges
#' 
#' A function for simulating the effect of incomplete sampling of the fossil
#' record.
#' 
#' @details 
#' This function implements a range of sampling models in continuous time. Be
#' default, sampling is simulated under the simplest model, where sampling
#' occurs as a Poisson process under a instantaneous sampling rate (r) which is
#' homogeneous through time and across lineages (Foote, 1997). Under this model,
#' the waiting times to sampling events are exponentially distributed, with an
#' average waiting time of 1/r. This useful property allows sampling to be
#' rapidly simulated for many taxa under this simple model in sampleRanges, by
#' repeatedly drawing waiting times between sampling events from an exponential
#' distribution. This is the model that is run when alpha, beta and
#' rTimeRatio are set to 1.
#' 
#' In addition to this simple model, sampleRanges also can consider a range of
#' additional models, including the hatP and incP options of Liow et al.
#' (2010). To describe the behavior of these models, users alter the default
#' values for alpha, beta and rTimeRatio. These parameters, and r, can either
#' be a single value which describes the behavior of the entire dataset or as a
#' vector, of same length as the number of taxa, which describes the per-taxon
#' value.  When any rTimeRatio, alpha or beta value is not equal to one, then
#' the sampling rate will vary across the duration of a taxon's temporal range.
#' In general, setting alpha and beta equal to a value above 2 produce "hat" or
#' bell-shaped curves where sampling rates peak at the midpoint of taxon
#' ranges, while setting them unequal will produce asymmetric bell curves
#' according to the beta function (Liow et al., 2010; Liow et al. set
#' alpha=beta=4). rTimeRatio is the ratio of the sampling rate of the
#' latest/most recent time divided by the earliest/oldest time.
#' 
#' The input r values will be interpreted differently based on whether one r
#' value or per-taxon values were used. If one value was input, then it is
#' assumed that r represent the grand mean r for the entire dataset for purposes
#' of time-varying r, such that if rTimeRatio is not equal to 1, taxa near the
#' end and start of the dataset will have very different per-taxon mean
#' sampling rate. If per-taxon values of r were input, then each r is consider
#' the per-taxon mean sampling rate. These will not be changed, but any
#' within-lineage variation is distributed so that the mean is still the input
#' per-taxon value. This also changes the interpretation of rTimeRatio, such
#' that when a single r value and rTimeRatio is given, it is assumed the ratio
#' describes the change in sampling rates from the start of the dataset to the
#' end, while if multiple values are given for either r or rTimeRatio will
#' instead see the value as describing the ratio at the first and last times of
#' each taxon. For the pure hat model, this interpretation of 'r' as a grand mean
#' sampling means that taxa will have a sampling rate of 2*r at the mid-peak of their
#' range, which will have considerable implications for taxonomic incompleteness.
#' 
#' The particular distinctions about these parameter values are important: all
#' models simulated in sampleRanges are structured to be effectively nested
#' inside a most general model with parameters r, alpha, beta and rTimeRatio.
#' 
#' Note that the modeling of sampling in this function is independent and
#' secondary of the actual simulation of the ranges, which are (generally)
#' produced by the models of simFossilRecord with argument \code{r}
#' (sampling rate) not set. Thus, 'hat-shaped range
#' distributions' are only contained within single morphotaxa; they do not
#' cross multiple morphotaxa in the case of anagenesis. Cryptic taxa each have
#' their own hat and do not share a single hat; by default the ranges of
#' cryptic taxa are merged to produce the range of a single observed
#' morphotaxon.
#' 
#' 'Hats' are constrained to start and end with a taxon's range, representing
#' the rise and fall of taxa in terms of abundance and geographic range (Liow
#' et al., 2010). However, for still-living taxa at the modern day, it is
#' unknown how much longer they may be alive (for memoryless Poisson models,
#' there is no age-dependent extinction). The treatment of these taxa with
#' regards to their 'hat' (i. e. the beta distribution) is controlled by randLivehat:
#' when ranLiveHat=FALSE, the beta distribution is fit so that the last appearance of
#' still-alive taxa at the modern day is treated as a last appearance for calculating the hat. When TRUE,
#' the default option, the still-alive taxa are considered to have gotten some
#' distance between 0 and 1 through the beta distribution, as of the modern
#' day. This point of progression is stochastically selected for each taxon by
#' pulling a number from a uniform distribution, and used for calculating the hat.
#' 
#' Because sampling rate varies over morphotaxon ranges under any of these more
#' complex models, sampling events cannot be quickly simulated as waiting times
#' pulled from an exponential distribution. Instead, the taxon durations are
#' discretized into a large number of small time intervals of length minInt
#' (see above; minInt should be small enough that only one sampling event could
#' feasibly happen per interval). The probability of an event occurring within
#' each interval is calculated and used to stochastically simulate sampling
#' events. For each interval, a number between 0 and 1 is randomly pulled from
#' a uniform distribution and to the per-interval sampling probability to test
#' if a sampling event occurred (if the random number is less than the
#' probability, a sampling event is recorded). In general, this method is
#' slower but otherwise comparable to the quicker waiting times method. See the
#' examples below for a small test of this.
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' If min.taxa is set to zero, the simulation may produce output in which no
#' taxa were ever sampled.
#' 
#' If modern.samp.prob is set to 1.0 (the default), then living taxa will
#' always be sampled at least at the present day (if there are any living
#' taxa). If the probability is less than 1, they will be sampled with that
#' probability at the modern day.
#' 
#' By default, this function will merge sampling events from morphologically
#' cryptic taxa, listing them as occurrences for the earliest member of that
#' group. To change this behavior, set merge.cryptic to FALSE.
#' 
#' Conditioning on sampling some minimum number of taxa may create strange
#' simulation results for some analyses, such as simulation analyses of
#' birth-death processes. Set min.taxa=0 to remove this conditioning.
#' 


#' @param taxad A two-column matrix of per-taxon ranges. The five-column matrix
#' output of \code{simFossilRecord}, post transformation with \code{fossilRecord2fossilTaxa}
#' can also be supplied, which will be common in simulation usages.

#' @param r Instantaneous average sampling rate per lineage time units; given
#' as a vector of length one or length equal to the number of taxa

#' @param alpha Alpha parameter of beta distribution; given as a vector of
#' length one or length equal to the number of taxa

#' @param beta Beta parameter of beta distribution; given as a vector of length
#' one or length equal to the number of taxa

#' @param rTimeRatio Ratio of most recent sampling rate over earliest sampling
#' rate; given as a vector of length one or length equal to the number of taxa

#' @param modern.samp.prob Probability of sampling living taxa at the present
#' day (time=0), see below.

#' @param min.taxa Minimum number of taxa sampled. The default is 2.

#' @param ranges.only If TRUE, gives taxon first and last occurrences only. If
#' FALSE, gives the time of all sampling events as a list.

#' @param minInt Minimum interval size used for simulating complex models

#' @param merge.cryptic If TRUE, sampling events for cryptic species will be
#' merged into one taxon.

#' @param randLiveHat If TRUE, taxa still alive at modern day have the
#' end-point of their 'hat' chosen from a uniform distribution.

#' @param alt.method If TRUE, use the alternative method of discretizing time
#' even if a simple model of sampling is being simulated.

#' @param plot If TRUE, plots the sampling models for each taxon against time.

#' @return If ranges.only is TRUE, then the output is a two-column per-taxon
#' matrix of first and last appearances in absolute time. NAs mean the taxon
#' was never sampled in the simulation.
#' 
#' If ranges.only is FALSE (the default), the output is a list, where each
#' element is a vector of sampling events the timing of sampling events, each
#' corresponding to a different taxon in the input. Elements that are NA are
#' unsampled taxa.
#' @author David W. Bapst
#' @seealso \code{\link{simFossilRecord}}, \code{\link{binTimeData}}
#' @references Foote, M. 1997 Estimating Taxonomic Durations and Preservation
#' Probability. \emph{Paleobiology} \bold{23}(3):278--300.
#' 
#' Liow, L. H., T. B. Quental, and C. R. Marshall. 2010 When Can Decreasing
#' Diversification Rates Be Detected with Molecular Phylogenies and the Fossil
#' Record? \emph{Systematic Biology} \bold{59}(6):646--659.
#' @examples
#' 
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' layout(1:2)
#' #let's see what the 'true' diversity curve looks like in this case
#' taxicDivCont(taxa)
#' #simulate a fossil record with imperfect sampling with sampleRanges
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #plot the diversity curve based on the sampled ranges
#' taxicDivCont(rangesCont)
#' #compare the true history to what we might observe!
#' 
#' #let's try more complicated models!
#' 
#' #a pull-to-the-recent model with x5 increase over time similar to Liow et al.'s  incP
#' layout(1:2)
#' rangesCont1 <- sampleRanges(taxa,r=0.5,rTimeRatio=5,plot=TRUE)
#' taxicDivCont(rangesCont1)
#' 
#' #a hat-shaped model
#' layout(1:2)
#' rangesCont1 <- sampleRanges(taxa,r=0.5,alpha=4,beta=4,plot=TRUE)
#' taxicDivCont(rangesCont1)
#' 
#' #a combination of these
#' layout(1:2)
#' rangesCont1 <- sampleRanges(taxa,r=0.5,alpha=4,beta=4,rTimeRatio=5,plot=TRUE)
#' taxicDivCont(rangesCont1)
#' 
#' #testing with cryptic speciation
#' layout(1)
#' recordCrypt<-simFossilRecord(p=0.1, q=0.1, prop.cryptic=0.5, nruns=1,
#'	nTotalTaxa=c(20,30), nExtant=0)
#' taxaCrypt<-fossilRecord2fossilTaxa(recordCrypt)
#' rangesCrypt <- sampleRanges(taxaCrypt,r=0.5)
#' taxicDivCont(rangesCrypt)
#' 
#' #an example of hat-shaped models (beta distributions) when there are live taxa
#' set.seed(444)
#' recordLive<-simFossilRecord(p=0.1, q=0.05, nruns=1,
#'	nTotalTaxa=c(5,100),nExtant=c(10,100))
#' taxaLive<-fossilRecord2fossilTaxa(recordLive)
#' #with end-points of live taxa at random points in the hat
#' rangesLive<-sampleRanges(taxaLive,r=0.1,alpha=4,beta=4,randLiveHat=TRUE,plot=TRUE)
#' #with all taxa end-points at end-point of hat
#' rangesLive<-sampleRanges(taxaLive,r=0.1,alpha=4,beta=4,randLiveHat=FALSE,plot=TRUE)
#' 
#' 
#' \donttest{
#' #simulate a model where sampling rate evolves under brownian motion
#' tree<-taxa2phylo(taxa,obs=taxa[,3])
#' sampRateBM <- rTraitCont(tree)
#' sampRateBM <- sampRateBM-min(sampRateBM)
#' layout(1:2)
#' rangesCont1 <- sampleRanges(taxa,r=sampRateBM,plot=TRUE)
#' taxicDivCont(rangesCont1)
#' 
#' #evolving sampling rate, hat model and pull of the recent
#' layout(1:2)
#' rangesCont1 <- sampleRanges(taxa,r=sampRateBM,alpha=4,beta=4,rTimeRatio=5,plot=TRUE)
#' taxicDivCont(rangesCont1)
#' layout(1)
#' 
#' #the simpler model is simulated by pulling waiting times from an exponential
#' #more complicated models are simulated by discretizing time into small intervals
#' #are these two methods comparable?
#' #let's look at the number of taxa sampled under both methods
#' summary(replicate(100,sum(!is.na(sampleRanges(taxa,r=0.5,alt.method=FALSE)[,1]))))
#' summary(replicate(100,sum(!is.na(sampleRanges(taxa,r=0.5,alt.method=TRUE)[,1]))))
#' #they look pretty similar!
#' }
#' 
#' @export sampleRanges
sampleRanges<-function(taxad,r,alpha=1,beta=1,rTimeRatio=1,modern.samp.prob=1,min.taxa=2,
	ranges.only=TRUE,minInt=0.01,merge.cryptic=TRUE,randLiveHat=TRUE,alt.method=FALSE,plot=FALSE){
	#sample ranges using a taxad matrix as input 
	#if (ranges.only=TRUE): outputs matrix of FADs/LADs, with NAs for unsampled taxa
	#if (ranges.only=FALSE): outputs per-species list with vectors of dates where that species was sampled
		#ranges and occurance are output on a BACKWORD-moving time-scale as expected for paleo data
	#if modern.samp.prob=1, then all still-living taxa (taxa at 0 for LAD) are ALWAYS last observed at zero
		#this approximates the fact that we think the present-day living biota is almost perfectly sampled
			#(well, relative to the modern)
	#if r is a vector, then it is considered to be a vector of per-species sampling rates
		#r is the mean sampling rate for a species
	#rTimeRatio is the ratio of the latest taxon-avg sampling rates by the earliest taxonavg sampling rates
		#i.e. the proportional increase over (a) the entire taxon's history or (b) the taxon duration (if specific)
		#either a single value for all taxa or taxon-specific values
	#all parameters can be given as single values or species-specific values
	#names<-paste("t",1:4,sep="");taxad<-cbind(c(250,230,210,200),c(240,215,205,0))
	#min.taxa=2;minInt=0.01;modern.samp.prob=0;plot=T;ranges.only=F;alt.method=F;randLiveHat=TRUE;merge.cryptic=TRUE
	#r<-c(0.2,0.1,0.3,0.4);alpha<-4;beta<-4;rTimeRatio<-2
	#r<-c(0,0.1,0.3,0.4);alpha<-beta<-rTimeRatio<-2
	if(ncol(taxad)==6){				#also allow it to accept taxad objects
		living<-taxad[,5]
		cryptic<-sapply(1:nrow(taxad),function(x) if(taxad[x,1]!=taxad[x,6]){which(taxad[,1]==taxad[x,6])}else{NA})
		timeData<-taxad[,3:4,drop=FALSE]
		names<-if(is.null(rownames(taxad))){paste("t",taxad[,1],sep="")}else{rownames(taxad)}
		}
	if(ncol(taxad)==2){			#assumes it has two matrices
		living<-rep(0,nrow(taxad))
		living[taxad[,2]==0]<-1
		cryptic<-rep(NA,nrow(taxad))
		timeData<-taxad
		names<-if(is.null(rownames(taxad))){paste("t",1:nrow(taxad),sep="")}else{rownames(taxad)}
		}
	if(nrow(taxad)<min.taxa){stop("min.taxa set higher than number of taxa in input")}	
	if(merge.cryptic==TRUE & nrow(taxad)<sum(!is.na(cryptic))){
		stop("min.taxa set higher than number of non-cryptic taxa in input")}
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Some dates in timeData <0 ?")}
	#check input parameters for inconsistencies
	if(any(r<0)){stop("some r values < 0")}
	if(any(alpha<=0) | any(beta<=0)){stop("some shape parameters (alpha, beta) <= 0")}
	if(length(alpha)==1){alpha<-rep(alpha,length(names))
		}else{if(length(alpha)!=length(names)){
		stop("Multiple alpha values input but not same length as number of taxa?")}}
	if(length(beta)==1){beta<-rep(beta,length(names))
		}else{if(length(beta)!=length(names)){
		stop("Multiple beta values input but not same length as number of taxa?")}}
	if(length(rTimeRatio)!=1){if(length(rTimeRatio)!=length(names)){
		stop("Multiple rTimeRatio values input but not same length as number of taxa?")}}
	midpoint<-(max(timeData)+min(timeData))/2
	taxa.midpoints<-((timeData[,1]+timeData[,2])/2)
	#need to calculate rTimeChange based on rTimeRatio
	r_start<-r/((rTimeRatio/2)+0.5)
	r_end<-2*r*(1-1/(rTimeRatio+1))
	if(length(r)==1 & length(rTimeRatio)==1){
		#thanks to Emily King and Brian Koch for helping me with this bit!
		if(max(timeData)==min(timeData)){
			denomTime<-1
			}else{
			denomTime<-(max(timeData)-min(timeData))
			}
		rTimeChange<-(r_end-r_start)/denomTime
		rTimeChange<-rep(rTimeChange,length(names))
	}else{	#if there's multiple r values, multiple ratios or both...
		dur1<-timeData[,1]-timeData[,2]
		rTimeChange<-(r_end-r_start)/dur1
		}
	if(length(r)==1){		#get per-taxon r.avg for each taxon, if not input already
		taxa.timeChange<-midpoint-taxa.midpoints
		r<-r+(rTimeChange*taxa.timeChange)
		}
	if(length(r)!=length(names)){stop("Multiple r values given but not same length as number of taxa?")}
	if(any(r<0)){stop("rTimeRatio so low as to cause some species-specific avg r < 0")}
	if(any(alpha!=1) | any(beta!=1) | any(rTimeChange!=0) | alt.method){		#get per species time vectors of r
		rangesTimes<-apply(timeData,1,function(x) seq(x[1],x[2],by=-minInt))
		#get the time change component per minInt of the taxon ranges
		#first get the different of each min int from the midpoints, then calculate the change r accordint to rTimeChange
		taxa.shiftTimes<-lapply(1:length(r),function(x) taxa.midpoints[x]-rangesTimes[[x]])
		taxa.rShiftTime<-lapply(1:length(r),function(x) rTimeChange[x]*taxa.shiftTimes[[x]])
		#if randLiveHat is TRUE, scale 
		#for each taxon decide on an end.hat; use runif if extant, 1 is extinct
		if(randLiveHat){end.hat<-ifelse(living==1,runif(length(living)),1)
			}else{end.hat<-rep(1,length(rangesTimes))}
		#now get the hat component, combine with the time component
		scaledTimes<-lapply(1:length(end.hat),function(x) seq(0,end.hat[x],length.out=length(rangesTimes[[x]])))
		rHat<-lapply(1:length(r),function(x) r[x]*dbeta(scaledTimes[[x]],alpha[x],beta[x]))
		rHatTime<-lapply(1:length(r),function(x) rHat[[x]]+taxa.rShiftTime[[x]])
		#transform so that the hats always are above zero but have correct means and correct timechange slope (so complicated!)
		#rHatTime<-lapply(rHatTime,function(x) (x-min(x))/(1-(min(x)/mean(x))))
		rHatTime1<-lapply(rHatTime,function(x) ifelse(x<0,0,x))
		rHatTime<-lapply(1:length(rHatTime),function(x) if(mean(rHatTime1[[x]])>0){
				mean(rHatTime[[x]])*rHatTime1[[x]]/mean(rHatTime1[[x]])
			}else{mean(rHatTime[[x]])*rHatTime1[[x]]}
			)
	}else{	
		rHatTime<-NULL
		}
	if(plot){if(!is.null(rHatTime)){
		plot(0:1,c(min(unlist(rHatTime)),max(unlist(rHatTime))),type="n",xlim=c(max(timeData),min(timeData)),	
			xlab="Time (Time-Units before Present)",ylab=c("Instant. Sampling Rate","(per lineage time-units)"))
		for(i in 1:length(r)){lines(rangesTimes[[i]],rHatTime[[i]],lwd=2)}
	}else{
		plot(c(max(timeData),min(timeData)),c(min(r),max(r)),type="n",xlim=c(max(timeData),min(timeData)),
			xlab="Time (Time-Units before Present)",ylab=c("Instant. Sampling Rate","(per lineage time-units)"))
		for(i in 1:length(r)){lines(timeData[i,],c(r[i],r[i]),lwd=2)}
		}}
	#plot(rangesTimes[[2]],rHatTime[[2]],type="l",lwd=3,xlim=c(max(rangesTimes[[2]]),
	#	min(rangesTimes[[2]])),xlab="time",ylab="Instantaneous Sampling Rate (per Lmy)")
	#(x<-cbind(r,sapply(rHatTime,mean)));plot(x);abline(0,1)	#means should be close to actual r
	#okay, now to actually do the sampling
	redo<-TRUE
	while(redo){
		samp_occ<-list()	#sampled occurances
		for(i in 1:nrow(timeData)){
			if(is.null(rHatTime)){
				if(r[i]>0){
					samps<-timeData[i,1]	#the time of speciation is the lower bound
					while(min(samps)>timeData[i,2]){	#keep sampling until you go past the time of extinction
						samps<-c(samps,min(samps)-rexp(1,rate=r[i]))}
					samps<-samps[-c(1,length(samps))]
				}else{samps<-numeric(0)}
			}else{
				samps<-rangesTimes[[i]][sapply(rHatTime[[i]],function(x) runif(1)<=(x*minInt))]
				}
			if(modern.samp.prob>0){if(living[i]==1){		#rework modern.sampling as a probability
				if(runif(1)<=modern.samp.prob){samps<-c(samps,0)}
				}}
			if(length(samps)>0){samp_occ[[i]]<-samps}else{samp_occ[[i]]<-NA}
			}
		redo<-min.taxa>sum(sapply(samp_occ,function(x) !any(is.na(x))))
		}
	names(samp_occ)<-names
	if(any(!is.na(cryptic)) & merge.cryptic){for(i in which(!is.na(cryptic))){
		samp_occ[[cryptic[i]]]<-c(samp_occ[[cryptic[i]]],samp_occ[[i]])
		samp_occ[[i]]<-NA
		}}
	if(ranges.only){
		ranges<-cbind(sapply(samp_occ,max),sapply(samp_occ,min))
		rownames(ranges)<-names;colnames(ranges)<-c("FAD","LAD")
		res<-ranges
	}else{
		res<-samp_occ
		}
	return(res)
	}
