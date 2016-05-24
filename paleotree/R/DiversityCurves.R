#' Diversity Curves
#' 
#' Functions to plot diversity curves based on taxic range data, in both
#' discrete and continuous time, and for phylogenies.
#' 
#' @details First, some background. Diversity curves are plots of species/taxon/lineage richness
#' over time for a particular group of organisms. For paleontological studies,
#' these are generally based on per-taxon range data while more recently in
#' evolutionary biology, molecular phylogenies have been used to calculate
#' lineage-through-time plots (LTTs). Neither of these approaches are without
#' their particular weaknesses; reconstructing the true history of biodiversity
#' is a difficult task no matter what data is available.
#' 
#' The diversity curves produced by these functions will always measure
#' diversity within binned time intervals (and plot them as rectangular bins).
#' For continuous-time data or phylogenies, one could decrease the int.length
#' used to get what is essentially an 'instantaneous' estimate of diversity.
#' This is warned against, however, as most historical diversity data will have
#' some time-averaging or uncertain temporal resolution and thus is probably
#' not finely-resolved enough to calculate instantaneous estimates of
#' diversity.
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' As diversity is counted within binned intervals, the true standing diversity
#' may be somewhat lower than the measured (observed) quantity, particularly if
#' intervals are longer than the mean duration of taxa is used. This will be an
#' issue with all diversity curve functions, but particularly the discrete-time
#' variant. For diversity data in particularly large discrete time intervals,
#' plotting this data in smaller bins which do not line up completely with the
#' original intervals will create a 'spiky' diversity curve, as these smaller
#' intersecting bins will have a large number of taxa which may have been
#' present in either of the neighboring intervals. This will give these small
#' bins an apparently high estimated standing diversity. This artifact is
#' avoided with the default setting split.int=TRUE, which will split any input
#' or calculated intervals so that they start and end at the boundaries of the
#' discrete-time range bins.
#' 
#' The timeList object should be a list composed of two matrices, the first
#' matrix giving by-interval start and end times (in absolute time), the second
#' matrix giving the by-taxon first and last appearances in the intervals
#' defined in the first matrix, numbered as the rows. Absolute time should be
#' decreasing, while the intervals should be numbered so that the number
#' increases with time. Taxa alive in the modern should be listed as last
#' occurring in a time interval that  begins at time 0 and ends at time 0.
#' See the documentation for the time-scaling  function 
#'\code{\link{bin_timePaleoPhy}} and the simulation function 
#' \code{\link{binTimeData}} for more information on formatting.
#'
#' Unlike some paleotree functions, such as  perCapitaRates, the intervals
#' can be overlapping or of unequal length. The diversity curve functions
#' deal with such issues by assuming taxa occur from the base of the interval
#' they are first found in until the end of the last interval they are occur
#' in. Taxa in wide-ranging intervals that contain many others will be treated
#' as occurring in all nested intervals. 
#' 
#' phyloDiv will resolve polytomies to be dichotomous nodes separated by
#' zero-length branches prior to calculating the diversity curve. There is no
#' option to alter this behavior, but it should not affect the use of the
#' function because the addition of the zero-length branches should produce an
#' identical diversity history as a polytomy. phyloDiv will also drop
#' zero-length terminal branches, as with the function dropZLB. This the
#' default behavior for the function but can be turned off by setting the
#' argument drop.zlb to FALSE.
#' 
#' 
#' @name DiversityCurves
#' @rdname DiversityCurves
#' @aliases taxicDivCont taxicDivDisc phyloDiv

#' @param timeData Two-column matrix giving the per-taxon first and last
#' appearances in absolute time. The simulated data tables output by \code{fossilRecord2fossilTaxa}
#' following simulation with \code{simFossilRecord} can also be supplied to \code{taxicDivCont}.

#' @param timeList A list composed of two matrices, giving interval start and end 
#' dates and taxon first and last occurrences within those intervals. See details.

#' @param tree A time-scaled phylogeny of class phylo.

#' @param int.length The length of intervals used to make the diversity curve.
#' Ignored if int.times is given.

#' @param int.times An optional two-column matrix of the interval start and end
#' times for calculating the diversity curve. If NULL, calculated internally.
#' If given, the argument split.int and int.length are ignored.

#' @param plot If TRUE, a diversity curve generated from the data is plotted.

#' @param plotLogRich If TRUE, taxic diversity is plotted on log scale.

#' @param drop.cryptic If TRUE, cryptic taxa are merged to form one taxon for
#' estimating taxon curves. Only works for objects from \code{simFossilRecord}
#' via \code{fossilRecord2fossilTaxa}.

#' @param drop.singletons If TRUE, taxa confined to a single interval will be
#' dropped prior to the diversity curve calculation. This is sometimes done if
#' single intervals have overly high diversities due to the 'monograph' effect
#' where more named taxa are known in certain intervals largely due to
#' taxonomic expert effort and not real changes in historical biotic diversity.

#' @param timelims Limits for the x (time) axis for diversity curve plots. Only
#' affects plotting. Given as either NULL (the default) or as a vector of
#' length two as for 'xlim' in the basic R function plot. Time axes will be plotted
#' \emph{exactly} to these values.

#' @param extant.adjust Amount of time to be added to extend start time for
#' (0,0) bins for extant taxa, so that the that 'time interval' doesn't appear 
#' to have an infinitely small width.

#' @param split.int For discrete time data, should calculated/input intervals
#' be split at discrete time interval boundaries? If FALSE, can create apparent
#' artifacts in calculating the diversity curve. See below.

#' @param drop.ZLB If true, zero-length terminal branches are dropped from the
#' input tree for phylogenetic datasets, before calculating standing diversity.

#' @return These functions will invisibly return a three-column matrix, where
#' the first two columns are interval start and end times and the third column
#' is the number of taxa/lineages counted in that interval.

#' @author David W. Bapst

#' @seealso \code{\link{multiDiv}}, \code{\link{timeSliceTree}},
#' \code{\link{binTimeData}}
#' 
#' There are several different functions for traditional LTT plots
#' (phylogenetic diversity curves), such as the function
#' ,\code{\link{ltt.plot}} in the package ape, the function \code{ltt} in the
#' package phytools, the function \code{plotLtt} in the package laser and the
#' function \code{LTT.average.root} in the package TreeSim.

#' @examples
#' 
#' #taxicDivDisc with the retiolinae dataset
#' data(retiolitinae)
#' taxicDivDisc(retioRanges)
#'
#' #simulation examples
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #let's see what the 'true' diversity curve looks like in this case
#' #plot the FADs and LADs with taxicDivCont()
#' taxicDivCont(taxa)
#' #simulate a fossil record with imperfect sampling with sampleRanges
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #plot the diversity curve based on the sampled ranges
#' layout(1:2)
#' taxicDivCont(rangesCont)
#' #Now let's use binTimeData to bin in intervals of 1 time unit
#' rangesDisc <- binTimeData(rangesCont,int.length=1)
#' #plot with taxicDivDisc
#' taxicDivDisc(rangesDisc)
#' #compare to the continuous time diversity curve
#' 
#' layout(1)
#' #Now let's make a tree using taxa2phylo
#' tree <- taxa2phylo(taxa,obs_time=rangesCont[,2])
#' phyloDiv(tree)
#' 
#' #a simple example with phyloDiv
#'   #using a tree from rtree in ape
#' set.seed(444)
#' tree <- rtree(100)
#' phyloDiv(tree)
#' 
#' #a neat example of using phyDiv with timeSliceTree 
#'  #to simulate doing molecular-phylogeny studies 
#'  #of diversification...in the past
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' taxicDivCont(taxa)
#' #that's the whole diversity curve
#' #with timeSliceTree we could look at the lineage accumulation curve 
#'  #we'd get of species sampled at a point in time
#' tree <- taxa2phylo(taxa)
#' #use timeSliceTree to make tree of relationships up until time=950 
#' tree950 <- timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=FALSE)
#' #use drop.extinct=T to only get the tree of lineages extant at time=950
#' tree950 <- timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=TRUE)
#' #now its an ultrametric tree with many fewer tips...
#' #lets plot the lineage accumulation plot on a log scale
#' phyloDiv(tree950,plotLogRich=TRUE)
#' 
#' #an example of a 'spiky' diversity curve and why split.int is a good thing
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' taxaDiv <- taxicDivCont(taxa)
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' rangesDisc <- binTimeData(rangesCont,int.length=10)
#' #now let's plot with taxicDivDisc() but with the intervals from taxaDiv
#'  #by default, split.int=TRUE
#' taxicDivDisc(rangesDisc,int.times=taxaDiv[,1:2],split.int=TRUE)
#' #look pretty
#' #now let's turn off split.int
#' taxicDivDisc(rangesDisc,int.times=taxaDiv[,1:2],split.int=FALSE)
#' #looks 'spiky'!
#' 
#' @export
taxicDivCont<-function(timeData,int.length=1,int.times=NULL,plot=TRUE,plotLogRich=FALSE,timelims=NULL,drop.cryptic=FALSE){
	#This function estimates diversity for bins from continuous-time range data
	#input is a per-species matrix of backwards-time FADs and LADs in 2 columns (FADs first)
		#assumes time is in millions of years
	#time interval starts and ends can be pre-input as a 2 column matrix
		#int.length is ignored in this case
	#output (if TRUE) is matrix of bin-start, bit-end, div
	tblen<-int.length
	if(ncol(timeData)==6){	#also allow it to accept taxad objects
		if(!drop.cryptic){
			timeData<-timeData[,3:4,drop=FALSE]
		}else{
			timeDataF<-sapply(unique(timeData[,6]),function(x) max(timeData[x==timeData[,6],3]))
			timeDataL<-sapply(unique(timeData[,6]),function(x) min(timeData[x==timeData[,6],4]))
			timeData<-cbind(timeDataF,timeDataL)
			}
		}	
	if(!inherits(timeData,"matrix")){
		if(inherits(timeData,"data.frame")){
			timeData<-as.matrix(timeData)
		}else{
			stop("timeData not of matrix or data.frame format")
			}
		}
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Some dates in timeData <0 ?")}
	FAD<-as.numeric(timeData[,1]);LAD<-as.numeric(timeData[,2])
	if(is.null(int.times)){
		midtimes<-seq(max(FAD)+2*tblen,min(LAD)-2*tblen,by=-tblen)
		midtimes<-midtimes[midtimes>0]
		int.start<-midtimes+(tblen/2)
		int.end<-midtimes-(tblen/2)
	}else{
		int.start<-int.times[,1];int.end<-int.times[,2]
		midtimes<-(int.start+int.end)/2
		}
	div<-sapply(1:length(midtimes),function(x) sum(FAD>=int.end[x])-sum(LAD>int.start[x]))
	if(plot){
		times1<-c(int.start,(int.end+((int.start-int.end)/10000)))
		div1<-c(div,div)[order(times1)]
		times1<-sort(times1)
		if(plotLogRich){
			plot(times1[div1>0],div1[div1>0],type="l",log="y",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xaxs=if(is.null(timelims)){"r"}else{"i"},
				xlab="Time (Before Present)",ylab="taxic Richness (Log Scale)")		
		}else{
			plot(times1,div1,type="l",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xaxs=if(is.null(timelims)){"r"}else{"i"},
				xlab="Time (Before Present)",ylab="Taxic Richness")
			}
		}
	res<-cbind(int.start,int.end,int.div=div)
	return(invisible(res))
	}

#' @rdname DiversityCurves
#' @export
taxicDivDisc<-function(timeList,int.times=NULL,drop.singletons=FALSE,plot=TRUE,plotLogRich=FALSE,timelims=NULL,
		extant.adjust=0.001,split.int=TRUE){
	#this function estimates diversity for binned intervals from discrete interval range data
	#input is a list with (1) interval times matrix and (2) species FOs and LOs
	#time interval starts and ends can be pre-input as a 2 column matrix
		#HOWEVER this could be pretty misleading!
		#standing richness may never be high as the apparent richness of some bins
	#output (if TRUE) is 3 col matrix of int-start, int-end, div
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
	#if(drop.extant){timeData[[2]][(timeData[[1]][timeData[[2]][,2],1]==0),1]<-NA}
	if(drop.singletons){timeData<-timeData[timeData[,1]!=timeData[,2],]}
	intMat[intMat[,1]==0,1]<-extant.adjust
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(!sapply(intMat,is.numeric))){stop("Some values in the interval times aren't numeric??")}
	if(any(apply(intMat,1,diff)>0)){stop("timeList[[1]] not in intervals in time relative to modern")}
	if(any(intMat[,2]<0)){stop("Some dates in timeList[[1]] <0 ?")}
	if(any(apply(timeData,1,diff)<0)){stop("timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Some dates in timeList[[2]] <0 ?")}
	Fint<-as.numeric(timeData[,1]);Lint<-as.numeric(timeData[,2])
	FAD<-intMat[Fint,1];LAD<-intMat[Lint,2]
	if(is.null(int.times)){
		avg_dur<-abs(mean(apply(timeList[[1]],1,diff)))
		int.bounds<-unique(c(intMat,max(intMat)+avg_dur,min(intMat)-avg_dur))		#add a little space at start and end
		int.bounds<-int.bounds[order(-int.bounds)]
		intMat<-cbind(int.bounds[-length(int.bounds)],int.bounds[-1])
		int.start<-intMat[,1];int.end<-intMat[,2]
		midtimes<-apply(intMat,1,mean)
	}else{
		if(split.int){	#if split.int, then any interval times given are split at discrete time intervals
			splinters<-sort(unique(c(intMat)))
			mustSplit<-apply(int.times,1,function(x) any(sapply(splinters,function(y) x[1]>y & x[2]<y)))
			if(any(mustSplit)){
				for(i in which(mustSplit)){
					splitter<-splinters[sapply(splinters,function(y) int.times[i,1]>y & int.times[i,2]<y)]
					for(j in splitter){		#in case there is more than one splitter
						int.times<-rbind(int.times,c(int.times[i,1],j),c(j,int.times[i,2]))
						}
					}
				int.times<-int.times[-which(mustSplit),]
				int.times<-int.times[order(-int.times[,1]),]
				}
			}
		int.start<-int.times[,1];int.end<-int.times[,2]
		midtimes<-(int.start+int.end)/2
		}
	div<-sapply(1:length(midtimes),function(x) sum(FAD>int.end[x])-sum(LAD>=int.start[x]))
	#div<-sapply(min(timeData):max(timeData),function(x) 	sum(FAD<=x & LAD>=x))
	if(plot){
		times1<-c(int.start,(int.end+((int.start-int.end)/10000)))
		div1<-c(div,div)[order(times1)]
		times1<-sort(times1)
		if(plotLogRich){
			plot(times1[div1>0],div1[div1>0],type="l",log="y",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xaxs=if(is.null(timelims)){"r"}else{"i"},
				xlab="Time (Before Present)",ylab="Taxic Richness (Log Scale)")		
		}else{
			plot(times1,div1,type="l",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xaxs=if(is.null(timelims)){"r"}else{"i"},
				xlab="Time (Before Present)",ylab="Taxic Richness")
			}
		}
	res<-cbind(int.start,int.end,int.div=div)
	return(invisible(res))
	}

#' @rdname DiversityCurves
#' @export
phyloDiv<-function(tree,int.length=0.1,int.times=NULL,plot=TRUE,plotLogRich=FALSE,
		drop.ZLB=TRUE,timelims=NULL){
	#function that computes a diversity curve from a tree file
		#aka lineage-through-time plot
	#root.time
		#ttree$root.time is used to place the diversity curve in time
		#if no root.time, then it is assumed latest tip is at 0 time (present day)
	#time interval starts and ends can be pre-input as a 2 column matrix
		#int.length is ignored in this case
	#this function will randomly resolve any tree it is given using multi2di()
		#this shouldn't affect anything to my knowledge
	#this function also automatically drops zero-length branches from the tree
		#this is advised for paleo-tree analyses of diversification
	#output (if TRUE) is 3 col matrix of bin-start, bit-end, div
	#plotLogRich just decides if the div plot if log-scale or not on the y axis
	#require(ape)
	ttree<-tree
	if(!inherits(ttree, "phylo")){
		stop("ttree is not of class phylo")
		}
	tblen<-int.length
	if(drop.ZLB){ttree<-dropZLB(ttree)}
	savetree<-ttree
	if(!is.binary.tree(ttree) | !is.rooted(tree)){ttree<-multi2di(ttree)}
	if(is.null(ttree$root.time)){
		ntime<-node.depth.edgelength(ttree)
		ntime<-max(ntime)-ntime
	}else{
		ntime<-node.depth.edgelength(ttree)
		ntime<-ttree$root.time-ntime
		ntime<-round(ntime,6)
		if(min(ntime)<0){stop("tree$root.time is less than total depth of tree!")}
		}
	if(is.null(int.times)){
		midtimes<-seq(max(ntime)+3*tblen,min(ntime)-2*tblen,by=-tblen)
		midtimes<-midtimes[midtimes>0]
		int.start<-midtimes+(tblen/2)
		int.end<-midtimes-(tblen/2)
	}else{
		int.start<-int.times[,1];int.end<-int.times[,2]
		midtimes<-(int.start+int.end)/2
		}
	LAD<-ntime[1:Ntip(ttree)]				#death
	FAD<-ntime[(Ntip(ttree)+1):length(ntime)]		#birth
	div<-sapply(1:length(midtimes),function(x) 1+sum(FAD>=int.end[x])-sum(LAD>int.start[x]))
	if(plot){
		times1<-c(int.start,(int.end+((int.start-int.end)/10000)))
		div1<-c(div,div)[order(times1)]
		times1<-sort(times1)
		layout(matrix(1:2,2,1))
		parOrig<-par(no.readonly=TRUE)
		par(mar=c(1,4,1,1))
		plot(ladderize(savetree),show.tip.label=FALSE)
		axisPhylo()    #anticipating that ape will recognize root.time soon
		par(mar=c(5,4,2,2))
		if(plotLogRich){
			plot(times1[div1>0],div1[div1>0],type="l",log="y",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xaxs=if(is.null(timelims)){"r"}else{"i"},
				xlab="Time (Before Present)",ylab="Lineage Richness (Log Scale)")		
		}else{
			plot(times1,div1,type="l",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xaxs=if(is.null(timelims)){"r"}else{"i"},
				ylim=c(0,max(div1)+1),
				xlab="Time (Before Present)",ylab="Lineage Richness")
			}
		par(parOrig)
		layout(1)
		}
	res<-cbind(int.start,int.end,int.div=div)
	return(invisible(res))
	}