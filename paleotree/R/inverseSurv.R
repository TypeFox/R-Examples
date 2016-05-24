#' Inverse Survivorship Models in the Fossil Record
#'
#' This function replicates the model-fitting procedure for forward and
#' reverse survivorship curve data, described by Michael Foote in a series
#' of papers (2001, 2003a, 2003b, 2005). These methods are discrete interval
#' taxon ranges, as are used in many other functions in paleotree 
#' (see function arguments). This function
#' can implement the continuous time, pulsed interval and mixed models described
#' in Foote (2003a and 2005).

#' @details
#' The design of this function to handle mixed continuous and discrete time models
#' means that parameters can mean very different things, dependent on how a
#' model is defined. Users should carefully evaluate their function arguments
#' and the discussion of parameter values as described in the Value section.

#' @inheritParams durationFreq

#' @param PA_n The probability of sampling a taxon after the last interval 
#' included in a survivorship study. Usually zero for extinct groups, 
#' although more logically has the value of 1 when there are still extant
#' taxa (i.e., if the last interval is the Holocene and the group is
#' still alive, the probability of sampling them later is probably 1...).
#' Should be a value of 0 to 1, NULL, or can be simply "fixed", the default option.
#' This default "fixed" option allows make_inverseSurv to decide the value
#' based on whether there is a modern interval (i.e. an interval that is 
#' \code{c(0,0)}) or not: if there is one, then \code{PA_n=1}, if not, 
#' then \code{PA_n=0}. If NULL, PA_n is treated as an additional free
#' parameter in the output model.

#' @param PB_1 The probability of sampling a taxon before the first interval 
#' included in a survivorship study. Should be a value of 0 to 1, or NULL. 
#' If NULL, PB_1 is treated as an additional free parameter in the output model.

#' @param p_cont If TRUE (the default), then origination is assumed to be a 
#' continuous time process with an instantaneous rate. If FALSE, the origination
#' is treated as a pulsed discrete-time process with a probability.

#' @param q_cont If TRUE (the default), then extinction is assumed to be a 
#' continuous time process with an instantaneous rate. If FALSE, the extinction
#' is treated as a pulsed discrete-time process with a probability.

#' @param Nb The number of taxa that enter an interval (b is for 'bottom'). This
#' is an arbitrary constant used to scale other values in these calculations and
#' can be safely set to 1.

#' @return
#' A function of class "paleotreeFunc", which takes a vector equal to the number
#' of parameters and returns the *negative* log likelihood (for use with optim and
#' similar optimizing functions, which attempt to minimize support values). See the
#' functions listed at \code{\link{modelMethods}} for manipulating and examining
#' such functions and \code{\link{constrainParPaleo}} for constraining parameters.
#'
#' The function output will take the largest number of parameters possible with
#' respect to groupings and time-intervals, which means the number of parameters
#' may number in the hundreds. Constraining the function for optimization
#' is recommended except when datasets are very large.
#'
#' Parameters in the output functions are named 'p', 'q' and 'r', which are
#' respectively the origination, extinction and sampling parameters. If the
#' respective arguments 'p_cont' and 'q_cont' are TRUE, then 'p' and 'q' will
#' represent the instantaneous per-capita origination and extinction rates
#' (in units of per lineage time-units). When one of these arguments is given as
#' FALSE, the respective parameter (p or q) will represent per-lineage-interval
#' rates. For p, this will be the per lineage-interval rate of a lineage producing
#' another lineage (which can exceed 1 because diversity can more than double) and
#' for q, this will be the per lineage-interval 'rate' of a lineage going extinct,
#' which cannot be observed to exceed 1 (because the proportion of diversity that
#' goes extinct cannot exceed 1). To obtain the per lineage-interval rates from a
#' set of per lineage-time unit rates, simply multiply the per lineage-time-unit
#' rate by the duration of an interval (or divide, to do the reverse; see Foote,
#' 2003 and 2005). 'r' is always the instantaneous per-capita sampling rate, in
#' units per lineage-time units. 
#'
#' If PA_n or PB_1 were given as NULL in the arguments, two additional parameters
#' will be added, named respectively 'PA_n' and 'PB_1', and listed separately for every
#' additional grouping. These are the probability of a taxon occurring before the first
#' interval in the dataset (PB_1) and the probability of a taxon occurring after
#' the last interval in a dataset (PA_n). Theses will be listed as 'PA_n.0' and 'PB_1.0'
#' to indicate that they are not related to any particular time-interval included
#' in the analysis, unlike the p, q, and r parameters (see below).
#'
#' Groupings follow the parameter names, separated by periods; by default, the
#' parameters will be placed in groups corresponding to the discrete intervals
#' in the input timeList, such that make_inverseSurv will create a function with
#' parameters 'p.1', 'q.1' and 'r.1' for interval 1; 'p.2', 'q.2' and 'r.2' for
#' interval 2 and so on. Additional groupings given by the user are listed after 
#' this first set (e.g. 'p.1.2.2').
#'
#' \subsection{Calculating The Results of an Inverse Survivorship Model}{
#' Because of the complicated grouping and time interval scheme, combined with
#' the probable preference of users to use constrained models rather that the
#' full models, it may be difficult to infer what the rates for particular
#' intervals and groups actually are in the final model, given the parameters
#' that were found in the final optimization.
#'
#' To account for this, the function output by \code{inverseSurv} also contains
#' an alternative mode which takes input rates and returns the final values along with
#' a rudimentary plotting function. This allows users to output per-interval and per-group
#' parameter estimates. To select these feature, the argument \code{altMode} must
#' be TRUE. This function will invisibly return the rate values for each
#' group and interval as a list of matrices, with each matrix composed of the
#' p, q and r rates for each interval, for a specific grouping.
#'
#' This plotting is extremely primitive, and most users will probably find the
#' invisibly returned rates to be of more interest. The function \code{layout} is
#' used to play plots for different groupings in sequence, and this may lead to
#' plots which are either hard to read or even cause errors (because of too many
#' groupings, producing impossible plots). To repress this, the argument \code{plotPar}
#' can be set to FALSE.
#'
#' This capability means the function has more arguments that just the
#' usual \code{par} argument that accepts the vector of parameters for running an
#' optimization. The first of these additional arguments, \code{altMode} enables
#' this alternative mode, instead of trying to estimate the negative log-likelihood
#' from the given parameters. The other arguments augment the calculation and plotting
#' of rates.
#'
#' To summarize, a function output by inverseSurv has the following arguments:
#' 
#' \describe{
#'  \item{par}{A vector of parameters, the same length as the number of parameters needed.
#' For plotting, can be obtained with optimization}

#'  \item{altMode}{If FALSE (the default) the function will work like ordinary model-fitting functions,
#' returning a negative log-likelihood value for the input parameter values in \code{par}. If TRUE,
#' however, the input parameters will instead be translated into the by-interval, by-group rates
#' used for calculating the log-likelihoods, plotted (if plotPar is TRUE) and these final
#' interval-specific rates will be returned invisibly as described above.}

#'  \item{plotPar}{If TRUE (the default) the calculated rates will be plotted, with each
#' grouping given a separate plot. This can be repressed by setting plotPar to FALSE. As the only
#' conceivable purpose for setting plotPar to FALSE is to get the calculated rates, these will not
#' be returned invisibly if plotPar is FALSE.}

#'  \item{ratesPerInt}{If FALSE, the default option, the rates plotted and returned will
#' be in units per lineage-time units, if those rates were being treated as rates for a
#' continuous-time process (i.e. p_cont=TRUE and q_cont=TRUE for p and q, respectively,
#' while r is always per lineage-time units). Otherwise, the respective rate will be in
#' units per lineage-interval. If ratesPerInt is TRUE instead, then \emph{all} rates, even
#' rates modelled as continuous-time process, will be returned as per lineage-interval rates,
#' even the sampling rate r.}

#'  \item{logRates}{If FALSE (the default) rates are plotted on a linear scale. If TRUE,
#' rates are plotted on a vertical log axis.}

#'  \item{jitter}{If TRUE (default) the sampling rate and extinction rate will be plotted slightly
#' ahead of the origination rate on the time axis, so the three rates can be easily differentiated. 
#' If false, this is repressed.}

#'  \item{legendPoisition}{The position of a legend indicating which line is
#' which of the three rates on the resulting plot. This
#' is given as the possible positions for argument 'x' of the function 
#' \code{\link{legend}}, and by default is "topleft", which will be generally
#' useful if origination and extinction rates are initially low. If 
#' legendPosition is NA, then a legend will not be plotted.}
#' }
#' }

#' @note This function is an entirely new rewrite of the methodology derived and presented by Foote in
#' his studies. Thus, whether it would give identical results cannot be assumed nor is it easy to test
#' given differences in the way data is handled between our coded functions. Furthermore, there may be
#' differences in the math due to mistakes in the derivations caught while this function was programmed.
#' I have tested the function by applying it to the same Sepkoski genus-level dataset that Foote used in
#' his 2003 and 2005 papers. Users can feel free to contact me for detailed figures from this analysis.
#' Overall, it seems my function captured the overall pattern of origination and sampling rates, at least
#' under a model where both origination and extinction are modeled as continuous-time processes. Extinction
#' showed considerably more variability relative to the published figures in Foote (2005). Additional
#' analyses are being run to identify the sources of this discrepancy, and the function is being released
#' here in paleotree on a trial basis, so that it can be more easily loaded onto remote servers. Users should be
#' thus forewarned of this essentially 'beta' status of this function.

#' @aliases make_inverseSurv invSurv

#' @seealso
#' This function extensively relies on \code{\link{footeValues}}.
#'
#' A similar format for likelihood models can be seen in \code{\link{durationFreq}}.
#'
#' Also see \code{\link{freqRat}}, \code{\link{sRate2sProb}},
#' \code{\link{qsRate2Comp}} \code{\link{sProb2sRate}} and \code{\link{qsProb2Comp}}.
#'
#' For translating between sampling probabilities and sampling rates, see
#' \code{\link{SamplingConv}}.

#' @author 
#' David W. Bapst, with some advice from Michael Foote.

#' @references
#' Foote, M. 2001. Inferring temporal patterns of preservation, origination, and 
#' extinction from taxonomic survivorship analysis. \emph{Paleobiology} 27(4):602-630.
#'
#' Foote, M. 2003a. Origination and Extinction through the Phanerozoic: A New
#' Approach. \emph{The Journal of Geology} 111(2):125-148.
#'
#' Foote, M. 2003b. Erratum: Origination and Extinction through the Phanerozoic:
#' a New Approach. \emph{The Journal of Geology} 111(6):752-753.
#'
#' Foote, M. 2005. Pulsed origination and extinction in the marine realm.
#' \emph{Paleobiology} 31(1):6-20.

#' @examples
#' \donttest{
#' 
#' # let's simulate some taxon ranges from an imperfectly sampled fossil record
#' set.seed(444)
#' record <- simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa <- fossilRecord2fossilTaxa(record)
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #bin the ranges into discrete time intervals
#' rangesDisc <- binTimeData(rangesCont,int.length=5)
#' 
#' #apply make_inverseSurv
#' likFun<-make_inverseSurv(rangesDisc)
#' #use constrainParPaleo to make the model time-homogenous
#'   	#match.all~match.all will match parameters so only 2 pars: p=q and r
#' constrFun<-constrainParPaleo(likFun,match.all~match.all)
#' results <- optim(parInit(constrFun), constrFun,
#'       lower=parLower(constrFun), upper=parUpper(constrFun),
#'       method="L-BFGS-B", control=list(maxit=1000000))
#' results
#'
#' #plot the results
#' constrFun(results$par, altMode=TRUE)
#'
#' #unconstrained function with ALL of 225 parameters!!!
#'     # this will take forever to converge, so it isn't run
#' optim(parInit(likFun),likFun,lower=parLower(likFun),upper=parUpper(likFun),
#'       method="L-BFGS-B",control=list(maxit=1000000))
#' }

#' @name inverseSurv
#' @rdname inverseSurv
#' @export
make_inverseSurv<-function(timeList,groups=NULL,p_cont=TRUE,q_cont=TRUE,
	PA_n="fixed",PB_1=0,Nb=1,drop.extant=TRUE){
		#
	#infer a seperate p, q, r for every interval
		#PAn can be free or constrained: if last interval is modern, 1, if last interval pre-modern, 0
		#PB1 can be free or constrained, in which case it is assumed to be 0 (no sampling before the
			#given time interval
	#groups should be given as matrix with integer codes identifying different groups
		#a different column for each group system
		#one row for each taxon in timeList, as entered (prior to dropModern)
	#dropModern, if TRUE drop all modern taxa
	#Nb is a nuisance parameter; all parameters scale to Nb, and usually set arbitrary to 1
	#examples
		#library(paleotree);set.seed(444)
		#record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
		#	nTotalTaxa=c(30,40), nExtant=0)
		#taxa<-fossilRecord2fossilTaxa(record)
		#rangesCont <- sampleRanges(taxa,r=0.5,,modern.samp.prob=1)
		#timeList <- binTimeData(rangesCont,int.length=1)
		#PA_n<-"fixed";p_cont=T;q_cont=F
		#groups<-cbind(sample(0:1,nrow(timeList[[2]]),replace=TRUE),
		#	sample(0:1,nrow(timeList[[2]]),replace=TRUE))
	#
		#groups=NULL;PA_n="fixed";PB_1=0;p_cont=TRUE;q_cont=TRUE;Nb=1;dropModern=TRUE
	#
	modernTest<-apply(timeList[[1]],1,function(x) all(x==0))
	if(any(modernTest)){	#if modern present
		if(sum(modernTest)>1){stop("More than one modern interval in timeList??!")}
		#modify the taxon occurrence matrix
		modInt<-which(modernTest)
		newInt<-which(apply(timeList[[1]],1,function(x) x[1]!=0 & x[2]==0))
		if(length(newInt)>1){stop("More than one interval stretching to the modern in timeList??!")}
		if(drop.extant){
			modDroppers<-timeList[[2]][,1]==modInt
			timeList[[2]]<-timeList[[2]][!modDroppers,]
			}
		#change all modInt references to the prior int in the taxon appearance matrix
		timeList[[2]]<-apply(timeList[[2]],2,sapply,function(x) if(x==modInt){newInt}else{x})
		#modify interval matrix
		timeList[[1]]<-timeList[[1]][-modInt,]
		if(PA_n=="fixed"){PA_n=1}
		}
	if(PA_n=="fixed"){PA_n=0}	#if no modern interval, PA_n when fixed is 0
	n<-nrow(timeList[[1]])	#number of intervals
	int.length<-(-apply(timeList[[1]],1,diff))
	#get parNames
	parNames<-c("p","q","r")
	parNames<-as.vector(sapply(parNames,function(x) paste(x,1:n,sep=".")))
	if(is.null(PB_1)){parNames<-c(parNames,c("PB_1.0"))}
	if(is.null(PA_n)){parNames<-c(parNames,c("PA_n.0"))}
	if(!is.null(groups)){
		if(drop.extant & any(modernTest)){groups<-groups[-modDroppers,]}
		if(nrow(timeList[[2]])!=nrow(groups)){
			stop(paste("number of rows in groups isn't equal to number of taxa in timeList",
				if(drop.extant){"after modern taxa are dropped"}))}
		for(i in 1:ncol(groups)){
			parNames<-as.vector(sapply(parNames,function(x) paste(x,unique(groups[,i]),sep=".")))
			}
		groupings<-unique(groups)
		}
	ngroup<-ifelse(is.null(groups),1,nrow(groupings))
	#break parnames into a character matrix
	breakNames<-t(sapply(parNames,function(x) unlist(strsplit(x,split=".",fixed=TRUE))))
	#set bounds
	lowerBound<-rep(0.001,length(parNames))
	upperBound<-rep(5,length(parNames))
	upperBound[breakNames[,1]=="PB_1"]<-1
	upperBound[breakNames[,1]=="PA_n"]<-1
	#01-06-14 p can exceed 1 because diversity can more than double
		# q cannot because each lineage can only go extinct once
	#if(!p_cont){upperBound[breakNames[,1]=="p"]<-1}
	if(!q_cont){upperBound[breakNames[,1]=="q"]<-1}
	parbounds<-list(lowerBound,upperBound)
	#
	#first build survivorship table
		#get Xij number of taxa that appeared in i and continued to interval j
		#need to get a different table for every grouping of taxa
	obs_X_list<-list()
	for(k in 1:ngroup){
		if(is.null(groups)){
			ranges<-timeList[[2]]
		}else{
			ranges<-timeList[[2]][apply(groups,1,function(x) all(x==groupings[k,])),]
			}
		#make survivorship table from selected ranges
		obs_X<-matrix(0,n,n)
		for(i in 1:n){for(j in 1:n){if(i<=j){
			obs_X[i,j]<-sum(ranges[,1]==i & ranges[,2]==j)
			}}}
		obs_X_list[k]<-list(obs_X)
		}
	#
	logL_invSurv<-function(par,altMode=FALSE,plotPar=TRUE,ratesPerInt=FALSE,
		logRates=FALSE,jitter=TRUE,legendPosition="topleft"){
		if(length(par)!=length(parNames)){stop("Number of input parameters is not equal to number of parnames")}
		if(!altMode){
			logLsum<-numeric()
			for(z in 1:ngroup){
				if(is.null(groups)){selector<-rep(TRUE,length(par))
					}else{selector<-apply(breakNames,1,function(x) x[-(1:2)]==groupings[z,])}
				parnew<-par[selector]
				breaknew<-breakNames[selector,]
				#get pars and order according to interval
				p<-(parnew[breaknew[,1]=="p"])[order(as.numeric(breaknew[breaknew[,1]=="p",2]))]
				q<-(parnew[breaknew[,1]=="q"])[order(as.numeric(breaknew[breaknew[,1]=="q",2]))]
				r<-(parnew[breaknew[,1]=="r"])[order(as.numeric(breaknew[breaknew[,1]=="r",2]))]
				#correct for int.length, make rates relative to interval length
					#this means the rates we input as par are (when they are continuous) rates per Ltu
					#taking product with interval length makes the rates per lineage interval
				r<-r*int.length
				if(p_cont){p<-p*int.length}
				if(q_cont){q<-q*int.length}
				#	}
				obs_X<-unlist(obs_X_list[[z]])
				#
				fR<-footeValues(p=p,q=q,r=r,PA_n=PA_n,PB_1=PB_1,p_cont=p_cont,q_cont=q_cont,Nb=Nb)
				#P_forw
				P_forw<-matrix(,n,n)
				for(j in 1:n){
					P_forw[j,j]<-fR$XFL[j]/(fR$XFt[j]+fR$XFL[j])
					}
				for(i in 1:n){for(j in 1:n){if(i<j){
					if((i+1)<=(j-1)){k<-(i+1):(j-1)}else{k<-NA}
					if(q_cont){
						P_forw[i,j]<-(1-P_forw[i,i])*ifelse(is.na(k[i]),1,exp(-sum(q[k])))*((1-exp(-q[j]))
							*fR$PD_bL[j]+exp(-q[j])*fR$PD_bt[j]*(1-fR$PA[j]))/fR$PA[i]
						}else{
							P_forw[i,j]<-(1-P_forw[i,i])*ifelse(is.na(k[i]),1,prod(1-q[k]))*(q[j]*
								fR$PD_bL[j]+(1-q[j])*fR$PD_bt[j]*(1-fR$PA[j]))/fR$PA[i]
						}
					}}}
				#P_inv
				P_inv<-matrix(,n,n)
				for(j in 1:n){
					P_inv[j,j]<-fR$XFL[i]/(fR$XbL[i]+fR$XFL[i])
					}
				for(i in 1:n){for(j in 1:n){if(i<j){
					if((i+1)<=(j-1)){k<-(i+1):(j-1)}else{k<-NA}
					if(q_cont){
						P_inv[i,j]<-(1-P_inv[j,j])*ifelse(is.na(k[1]),1,exp(-sum(p[k])))*((1-exp(-p[i]))
							*fR$PD_Ft[i]+exp(-p[i])*fR$PD_bt[i]*(1-fR$PB[i]))/fR$PB[j]
					}else{
						P_inv[i,j]<-(1-P_inv[j,j])*ifelse(is.na(k[1]),1,prod(1/(1+p[k])))*(p[i]/(1+p[i])
							*fR$PD_Ft[i]+1/(1+p[i])*fR$PD_bt[i]*(1-fR$PB[i]))/fR$PB[j]
						}
					}}}
				#get log likelihood
				logL<-matrix(NA,n,n)
				for(i in 1:n){for(j in 1:n){if(obs_X[i,j]>0){
					logL[i,j]<-obs_X[i,j]*log(P_forw[i,j])+obs_X[i,j]*log(P_inv[i,j])
					#print(paste(logL[i,j],i,j))
					}}}
				logLsum[z]<-sum(logL,na.rm=TRUE)
				}
				logLsum<-sum(logLsum)
			return(unname(-logLsum))	#negate so optim works correctly
		}else{
			if(plotPar){
				if(ngroup<6){
					layout(1:ngroup)	#first do layout
				}else{
					message("Number of groups greater than 5, not plotting the resulting rates.")
					plotPar<-FALSE
					}
				}
			#
			results<-list()
			for(z in 1:ngroup){
				if(is.null(groups)){
					selector<-rep(TRUE,length(par))
				}else{
					selector<-apply(breakNames,1,function(x) x[-(1:2)]==groupings[z,])}
				parnew<-par[selector]
				breaknew<-breakNames[selector,]
				#get pars and order according to interval
				p<-(parnew[breaknew[,1]=="p"])[order(as.numeric(breaknew[breaknew[,1]=="p",2]))]
				q<-(parnew[breaknew[,1]=="q"])[order(as.numeric(breaknew[breaknew[,1]=="q",2]))]
				r<-(parnew[breaknew[,1]=="r"])[order(as.numeric(breaknew[breaknew[,1]=="r",2]))]
				if(ratesPerInt){ 
					#correct for int.length, make rates per interval, not per time
					r<-r*int.length
					if(p_cont){p<-p*int.length}
					if(q_cont){q<-q*int.length}
					}
				rates<-cbind(p,q,r)
				rownames(rates)<-sapply(strsplit(rownames(rates),"p"),function(x) x[2])
				results[[z]]<-rates
				if(plotPar){
					#now plot
					int.start<-timeList[[1]][,1];int.end<-timeList[[1]][,2]
					times1<-c(int.start,(int.end+((int.start-int.end)/100)))
					#add a jigger so rates don't overlap
					jigger<-min(diff(times1))/10
					p1<-c(p,p)[order(times1)]
					q1<-c(q,q)[order(times1)]
					r1<-c(r,r)[order(times1)]
					times1<-sort(times1)
					if(logRates){
						p1[!(p1>0)]<-NA
						q1[!(q1>0)]<-NA
						r1[!(r1>0)]<-NA
						ylims<-c(min(c(p1,q1,r1),na.rm=TRUE),(max(c(p1,q1,r1),na.rm=TRUE))*1.5)
						plot(times1,p1,type="l",log="y",col=4,lwd=2,lty=5,
						xlim=c(max(times1),max(0,min(times1))),ylim=ylims,
						xlab="Time (Before Present)",ylab="Instantaneous Per-Capita Rate (per Ltu)")
					}else{
						ylims<-c(min(c(p1,q1,r1),na.rm=TRUE),(max(c(p1,q1,r1),na.rm=TRUE))*1.5)
						plot(times1,p1,type="l",col=4,lwd=2,lty=5,
						xlim=c(max(times1),max(0,min(times1))),ylim=ylims,
						xlab="Time (Before Present)",ylab="Instantaneous Per-Capita Rate (per Ltu)")
						}
					if(jitter){
						lines(times1+jigger,q1,col=2,lwd=2,lty=2)
						lines(times1-jigger,r1,col=3,lwd=2,lty=6)
					}else{
						lines(times1,q1,col=2,lwd=2,lty=2)
						lines(times1,r1,col=3,lwd=2,lty=6)
						}
					if(!is.na(legendPosition)){
						legend(x=legendPosition,
							legend=c("Origination","Extinction","Sampling"),
							lty=c(5,2,6),lwd=2,col=c(4,2,3))
						}
					}
				}
			if(ngroup==1){results<-results[[1]]}
			#return the table of broken recalculated rates
			#if(plotPar){
			#return(invisible(results))
			#}else{
			return(results)
			#}
			}
		}
	#make into a paleoFunc
	logL_invSurv<-make_paleotreeFunc(logL_invSurv,parNames,parbounds)
	return(logL_invSurv)
	}
	