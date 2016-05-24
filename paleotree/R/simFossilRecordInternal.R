
#functions that will live only in simFossilRecord's namespace
	#so not to crowd paleotree namespace

makeParFunct<-function(par,isBranchRate){
	#things to look for
		# N = number of extant taxa
		# T = time
		# D = duration
		# P = branching rate
	acceptedArg<-if(isBranchRate){c('N','T','D')
		}else{c('N','T','D','P')}
	if(is.numeric(par)){par<-as.character(par)}
	#set default timeDepAttr to FALSE
	timeDepAttr<-FALSE
	if(is.numeric(type.convert(par,as.is=TRUE))){
		res<-type.convert(par,as.is=TRUE)
		#convert any 'Inf' rates to 0: this event type cannot occur !
		if(is.infinite(res)){res<-0}
		#return error if any rate is negative
		if(res<0){
			stop("input rates must be at least 0 if input can be coerced to type numeric")}
		#now convert formula expression to function
		parFunct<-if(isBranchRate){
			function(N,T,D){}
		}else{
			function(N,T,D,P){}
			}
		body(parFunct)<-res
		res<-parFunct
	}else{
		#first convert par to a formula
		#check arguments match only accepted list
		args<-all.vars(as.formula(paste0('XXdave~',par)))
		args<-args[!(args=='XXdave')]
		if(!all(sapply(args,function(x) any(x==acceptedArg)))){
			if(isBranchRate){
				stop(paste0('Incorrect parameterization of branching rate formula, \n',
					'Only N and T allowed as variables'))
			}else{
				stop(paste0('Incorrect parameterization of a non-branching rate formula, \n',
					'Only P, N and T allowed as variables'))
				}
			}
		# define a logical object to be output as an attr that signifies
			# if a rate had "T" or "D" args	
		if(any(sapply(args,function(x) any(x==c("T","D"))))){
			timeDepAttr<-TRUE
			}
		#now convert formula expression to function
		parFunct<-if(isBranchRate){
			function(N,T,D){}
		}else{
			function(N,T,D,P){}
			}
		body(parFunct)<-parse(text=par)
		res<-parFunct
		}
	# add an attr saying whether or not its a time-dependent rate function
	attr(res,"timeDep")<-timeDepAttr
	return(res)
	}

#RANDOM EXAMPLES TO TEST makeParFunct with...
#makeParFunct(0.1,isBranchRate=TRUE)
#
#z<-0.1
#makeParFunct(z^2,isBranchRate=TRUE)
#
#makeParFunct('P-0.1*N',isBranchRate=TRUE)
#
#makeExtRate<-makeParFunct('P-0.1*N',isBranchRate=FALSE)
#makeExtRate(P=0.1,N=10,T=100)
#
#makeParFunct('0.1+T*0.2-0.1^N',isBranchRate=FALSE)

# get rate vector

getRateMatrix<-function(taxa,timePassed,taxaDurations,
		getBranchRate,getExtRate,getSampRate,getAnagRate,
		prop.cryptic,prop.bifurc,negRatesAsZero){
	#
	#get some basic summary statistics first
	#vector of which taxa are still alive
	whichExtant<-whichLive(taxa)
	#standing number of extant lineages 
	# (i.e. number of lineages that stuff can happen to)
	nLive<-length(whichExtant)	
	#
	###########################################################
	# calculate rates (which may be time or diversity dependent)
		#use rate-getting functions from above
	#
	# set up mega vector for all taxa
	rateMat<-matrix(,nLive,6)
	colnames(rateMat)<-c('budd','bifurc','anag','crypt','ext','samp')
	attr(rateMat,"whichExtant")<-whichExtant
	#get rates for each taxon
	for(i in 1:nLive){
		taxonDur<-taxaDurations[i]
		#get the new branching rate, extinction rate, sampling rate, anagenesis rate
		branchRate<-getBranchRate(N=nLive, T=timePassed, D=taxonDur)
		extRate<-getExtRate(N=nLive, T=timePassed, D=taxonDur, P=branchRate)
		sampRate<-getSampRate(N=nLive, T=timePassed, D=taxonDur, P=branchRate)
		anagRate<-getAnagRate(N=nLive, T=timePassed,  D=taxonDur, P=branchRate)
		##
		# now deal with proportional types of branching
		#get cryptic, budding and bifurcation components
		crypticRate<-branchRate*(prop.cryptic)
		#rate of morph differentiation per branching event
		morphRate<-branchRate*(1-prop.cryptic)
		buddRate<-morphRate*(1-prop.bifurc)
		bifurcRate<-morphRate*(prop.bifurc)		
		#
		#get probabilities of event types into rateVector
		rateVector<-c(buddRate,bifurcRate,anagRate,crypticRate,extRate,sampRate)	
		names(rateVector)<-c('budd','bifurc','anag','crypt','ext','samp')
		#check rates, make sure none are less than zero
		if(any(rateVector<0)){
			if(negRatesAsZero){
				rateVector[rateVector<0]<-0
			}else{
				stop(paste0(names(which(rateVector<0)),'rate calculated less than zero'))
				}
			}
		rateMat[i,]<-rateVector
		}
	# check that not *all* rates are 0
	if(!(sum(rateMat)>0)){
		stop("Simulation found scenario in which all rates at some time point are zero (?!)")
		}	
	#
	return(rateMat)
	}

# example of randomly sampling a matrix with weighted probs
	# and getting row/col indices
# m<-matrix(runif(6*6),ncol=6,nrow=6)
# m<-m/sum(m)
# m[1,2]<-1
# m<-m/sum(m)
# arrayInd(sample(length(m),1,prob=m),dim(m)) 
	

#internal functions for branching/extinction/anagenesis processes	

initiateTaxa<-function(startTaxa,time){
	newTaxa<-lapply(1:startTaxa,function(x) 
		newTaxon(newID=x,ancID=NA,time=time,looksLike=x)
		)
	return(newTaxa)
	}
	
newTaxon<-function(newID,ancID,time,looksLike){
	#creates an entirely new just-originated taxon
	#store taxa as a list structure
		# $taxa.data, exactly like 'taxa' output from fossilRecord2fossilTaxa
	taxaData<-c(newID,ancID,time,NA,1,looksLike)
	names(taxaData)<- c('taxon.id','ancestor.id','orig.time','ext.time','still.alive','looks.like')
	# $sampling.times = times of sampling events for this taxon
		#thus can come up with quick/simple ways of evaluating run conditions
		# e.g. evaluate number of sampled taxa by sum(length($sampling.times)>0)
	taxon<-list(taxa.data=taxaData, sampling.times=numeric())
	return(taxon)
	}

origination<-function(taxa,ancID,time,looksLike=NULL){
	#adds a new taxon via branching or anagenesis
	newID<-max(sapply(taxa,function(x) x[[1]][1]))+1
	if(is.null(looksLike)){
		looksLike<-newID
	}else{
		#looksLike needs to be looksLike of parent
		looksLike<-taxa[[looksLike]][[1]][6]
		}
	newTaxonData<-newTaxon(newID=newID,ancID=ancID,time=time,looksLike=looksLike)	
	newTaxa<-c(taxa,list(newTaxonData))
	return(newTaxa)
	}

termination<-function(taxa,target,time){
	#ends an existing taxon (from extinction or pseudoextinction)
	whichTarget<-which(sapply(taxa,function(x) x[[1]][1]==target))
	if(length(whichTarget)!=1){stop('taxon IDs repeated??')}
	taxa[[whichTarget]][[1]][4]<-time
	taxa[[whichTarget]][[1]][5]<-0
	return(taxa)
	}

buddingEvent<-function(taxa,parent,time){
	taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
	return(taxa)
	}

crypticEvent<-function(taxa,parent,time){
	taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=parent)
	return(taxa)
	}

anagEvent<-function(taxa,parent,time){
	taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
	taxa<-termination(taxa=taxa,target=parent,time=time)
	return(taxa)
	}	
		
bifurcEvent<-function(taxa,parent,time){
	taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
	taxa<-origination(taxa=taxa,ancID=parent,time=time,looksLike=NULL)
	taxa<-termination(taxa=taxa,target=parent,time=time)
	return(taxa)
	}	
	
extEvent<-function(taxa,target,time){
	taxa<-termination(taxa=taxa,target=target,time=time)
	return(taxa)		
	}
	
sampEvent<-function(taxa,target,time){
	whichTarget<-which(sapply(taxa,function(x) x[[1]][1]==target))
	taxa[[whichTarget]][[2]]<-c(taxa[[whichTarget]][[2]],time)
	return(taxa)		
	}
		
eventOccurs<-function(taxa,target,type,time){
	#possible types : 'budd','bifurc','anag','crypt','ext','samp', 'wait'
	if(type=="budd"){
		taxa<-buddingEvent(taxa=taxa,parent=target,time=time)
		}
	if(type=="bifurc"){
		taxa<-bifurcEvent(taxa=taxa,parent=target,time=time)		
		}
	if(type=="anag"){
		taxa<-anagEvent(taxa=taxa,parent=target,time=time)
		}
	if(type=="crypt"){
		taxa<-crypticEvent(taxa=taxa,parent=target,time=time)
		}
	if(type=="ext"){
		taxa<-extEvent(taxa=taxa,target=target,time=time)		
		}	
	if(type=="samp"){
		taxa<-sampEvent(taxa=taxa,target=target,time=time)
		}
	#check
	if(any(sapply(taxa,length)!=2)){
		stop("taxa object contains taxa with not 2 elements?")}
	return(taxa)
	}

# functions for identifying live/sampled taxa and checking simulation runs

whichLive<-function(taxa){
	res<-which(sapply(taxa,function(x) x[[1]][5]==1))
	res2<-which(sapply(taxa,function(x) is.na(x[[1]][4]) | identical(unname(x[[1]][4]),0)))
	if(!identical(unname(res),unname(res2))){
		#browser()
		stop("Disagreement on which taxa are extant")}
	return(res)
	}

whichSampled<-function(taxa){
	res<-which(sapply(taxa,function(x) length(x[[2]])>0))
	return(res)
	}

getTaxonDurations<-function(taxa,currentTime){
	areExtant<-whichLive(taxa)
	durations<-sapply(taxa[areExtant],function(x) x[[1]][3]-currentTime)
	if(any(durations<0)){
		stop("negative durations calculated??")
		}
	return(durations)
	}
	
getRunVitals<-function(taxa,count.cryptic){
	#NOTE need to change vital measurement dependent on count.cryptic or not
	if(!is.list(taxa)){stop("handed getRunVitals a taxa object that isn't a list??")}
	whichExtant<-whichLive(taxa)
	whichSamp<-whichSampled(taxa)
	if(count.cryptic){
		nTaxa<-length(taxa)	#total number of taxa
		nLive<-length(whichExtant)
		nSampled<-length(whichSamp)	#total number of sampled taxa
	}else{
		looksLike<-sapply(taxa,function(x) x[[1]][6])
		#count number of unique taxa based on looksLike
		nTaxa<-length(unique(looksLike))
		#count number of unique extant taxa
		nLive<-length(unique(looksLike[whichExtant]))
		#count number of unique sampled taxa
		nSampled<-length(unique(looksLike[whichSamp]))
		}
	vitals<-c(nTotalTaxa=nTaxa,nExtant=nLive,nSamp=nSampled)
	return(vitals)
	}

testContinue<-function(vitals,timePassed,runConditions){
	#(1) continue = TRUE until max totalTime, max nTotalTaxa, nSamp or total extinction
		# none of these can ever REVERSE
	#
	#time passed
	#timePassed<-runConditions$totalTime[2]-currentTime			
	#
	# test run conditions
	totalExtinction<-vitals[2]==0
	tooMuchTime<-timePassed>runConditions$totalTime[2]
	tooManyTaxa<-vitals[1]>runConditions$nTotalTaxa[2]
	tooManySamp<-vitals[3]>runConditions$nSamp[2]
	runStop<-totalExtinction | tooMuchTime | tooManyTaxa | tooManySamp			
	continue<-unname(!runStop)
	return(continue)
	}

contiguousIntegerSeq<-function(vector){
	if(!is.vector(vector)){
		stop("contiguousIntegerSeq handed a non-vector?")}
	if(length(vector)<2){
		stop("contiguousIntegerSeq handed length 1 vector?")}
	#		
	vector<-as.integer(vector)
	#because unbelievably base R has no simple function for
		#pulling contiguous sequences of integers from 
	starts<-sapply(2:length(vector),function(i){
		(vector[i]-vector[i-1])>1
		})
	#browser()
	starts<-c(TRUE,starts)
	starts<-vector[starts]
	ends<-sapply(1:(length(vector)-1),function(i){
		(vector[i+1]-vector[i])>1
		})
	ends<-c(ends,TRUE)
	ends<-vector[ends]
	seqMat<-cbind(starts,ends)
	#
	# checks
	if(!is.matrix(seqMat)){stop("seqMat isn't a matrix")}
	if(ncol(seqMat)!=2){stop("seqMat doesn't have 2 columns")}
	return(seqMat)
	}

insertRow<-function(table,row,rownum){
	#because unbelievably base R has no simple function for inserting a row
	# insert new immediately at this number, shifts row currently at that location DOWN
	table<-rbind(table[1:(rownum-1),],row,table[-(1:(rownum-1)),])
	return(table)
	}
	
worthCheckingVitalsRecord<-function(vitalsRecord,runConditions){
	#
	# check that labels for vitalsRecord and runConditions match
	labMatch<-colnames(vitalsRecord)[2:4]==names(runConditions)[2:4]
	if(!all(labMatch)){
		stop("runConditions and vitalsRecord objects are mislabeled/misordered")
		}
	#
	lastVitals<-vitalsRecord[nrow(vitalsRecord),]
	reachMin<-sapply(c(1,2,4),function(i){
		lastVitals[i]>=runConditions[[i]][1]
		})
	worthChecking<-all(reachMin)
	#
	return(worthChecking)
	}

testVitalsRecord<-function(vitalsRecord,runConditions,tolerance){
	#
	# check that labels for vitalsRecord and runConditions match
	labMatch<-colnames(vitalsRecord)[2:4]==names(runConditions)[2:4]
	if(!all(labMatch)){
		stop("runConditions and vitalsRecord objects are mislabeled/misordered")
		}
	#
	#first INSERT FAKE EVENTS INTO VITALS MAT
		# FOR MIN TIME AND MAX TIME
	#
	# for min time
	if(vitalsRecord[1,1]<runConditions$totalTime[1] 
		& vitalsRecord[nrow(vitalsRecord),1]>runConditions$totalTime[1]
		& all(vitalsRecord[,1]!=runConditions$totalTime[1])){
		#
		#what row to insert at
		whereInsert<-which(vitalsRecord[,1]>runConditions$totalTime[1])[1]
		newRow<-c(runConditions$totalTime[1],vitalsRecord[whereInsert-1,-1])
		#vitalsRecord<-insertRow(table=vitalsRecord,row=newRow,rownum=whereInsert-1)
		vitalsRecord<-rbind(vitalsRecord,newRow)
		}
	vitalsRecord<-vitalsRecord[order(vitalsRecord[,1]),]
	#
	# for max time
	if(vitalsRecord[1,1]<runConditions$totalTime[2] 
		& vitalsRecord[nrow(vitalsRecord),1]>runConditions$totalTime[2]
		& all(vitalsRecord[,1]!=runConditions$totalTime[2])){
		#
		#what row to insert at
		whereInsert<-rev(which(vitalsRecord[,1]<runConditions$totalTime[2]))[1]
		newRow<-c(runConditions$totalTime[2],vitalsRecord[whereInsert,-1])
		#vitalsRecord<-insertRow(table=vitalsRecord,row=newRow,rownum=whereInsert+1)
		vitalsRecord<-rbind(vitalsRecord,newRow)
		}
	vitalsRecord<-vitalsRecord[order(vitalsRecord[,1]),]
	#
	################################################################
	#
	# NOW need to essentially duplicate EVERY ROW with time-stamp of row immediately after it
	newVitalsRecord<-cbind(vitalsRecord[2:nrow(vitalsRecord),1,drop=FALSE],
		vitalsRecord[1:(nrow(vitalsRecord)-1),2:4,drop=FALSE])
	pastIncrement<-diff(vitalsRecord[,1])
	if(any(pastIncrement>0)){
		pastIncrement<-min(pastIncrement[pastIncrement>0])/1000
		pastIncrement<-min(c(tolerance,pastIncrement))
	}else{
		pastIncrement<-tolerance
		}
	newTimes<-newVitalsRecord[,1]-pastIncrement
	newTimes<-ifelse(newTimes>0,newTimes,0)
	newVitalsRecord[,1]<-newTimes
	vitalsRecord<-rbind(vitalsRecord,newVitalsRecord)
	vitalsRecord<-vitalsRecord[order(vitalsRecord[,1]),]
	#
	###########################################################################
	# identify all rows where timePassed, nTaxa, nExtant and nSamp are good
	#
	okayVitalsMat<-sapply(1:4,function(i){
		var<-vitalsRecord[,i]
		varRange<-runConditions[[i]]
			var>=varRange[1] & var<=varRange[2]
			})
	okayVitals<-apply(okayVitalsMat,1,all)				
	#
	#########################################################
	# Now test if there are any, if so, sequence
	#
	if(any(okayVitals)){
		# need to build a matrix of the paired-date sequences
		whichOkay<-which(okayVitals)
		if(length(whichOkay)>1){
			seqVitals<-contiguousIntegerSeq(whichOkay)
		}else{
			seqVitals<-matrix(whichOkay,1,2)
			}
		#replaced with the passedTime dates
		seqVitals<-apply(seqVitals,2,function(x) vitalsRecord[x,1])
		if(is.vector(seqVitals)){
			seqVitals<-matrix(seqVitals,,2)
			}
		#
		# checks
		if(!is.matrix(seqVitals)){stop("seqVitals isn't a matrix")}
		if(ncol(seqVitals)!=2){stop("seqVitals doesn't have 2 columns")}
	}else{
		seqVitals<-NA
		}
	#
	#check to make sure it makes sense
	if(length(seqVitals)==0){stop("seqVitals constructed incorrectly")}
	#
	return(seqVitals)
	}

sampleSeqVitals<-function(seqVitals){
	cumSumSeq<-cumsum(apply(seqVitals,1,diff))
	totalSum<-rev(cumSumSeq)[1]
	if(totalSum>0){
		placedDate<-runif(1)*totalSum
		findRow<-which(cumSumSeq>=placedDate)[1]
		earlierRowCumSum<-ifelse(findRow==1,0,cumSumSeq[findRow-1])
		date<-seqVitals[findRow,1]+placedDate-earlierRowCumSum
	}else{
		# no probability density to sample
		# randomly pick a row
		if(length(seqVitals[,1])>1){
			date<-sample(seqVitals[,1],1)
		}else{
			date<-seqVitals[,1]
			}
		}					
	return(date)
	}

getTaxaNames<-function(taxa){	
	#name each normal taxon as t + ID 
		#cryptic taxa are cryptic id + . taxon number within that complex
	cryptIDs<-sapply(taxa,function(x) x[[1]][6])
	taxonIDs<-sapply(taxa,function(x) x[[1]][1])
	whichCrypt<-sapply(cryptIDs,function(x) sum(x==cryptIDs)>1)
	newIDs<-sapply(1:length(taxonIDs),function(x){
		if(whichCrypt[x]){
			#find what number this taxon is, within the cryptic complex
			matchCrypt<-cryptIDs==cryptIDs[x]
			nCrypt<-which(taxonIDs[matchCrypt]==taxonIDs[x])-1
			res<-paste0(cryptIDs[x],".",nCrypt)
		}else{
			res<-taxonIDs[x]
			}
		return(res)
		})
	newIDs<-paste0("t",newIDs)
	return(newIDs)
	}

testFinal<-function(taxa,timePassed,runConditions,count.cryptic){
	# need to adjust sampling for modern.sampling
	whichExtant<-whichLive(taxa)
	if(length(whichExtant)>0){
		# find extant time
		extantTime<-taxa[[whichExtant[1]]][[1]][4]
		for(i in 1:length(taxa)){
			# remove all sampling events at zero
			taxa[[i]][[2]]<-taxa[[i]][[2]][taxa[[i]][[2]]!=extantTime]
			}
		}
	# test that the produced taxa object actually passed the runConditions
	finalVitals<-getRunVitals(taxa=taxa,count.cryptic=count.cryptic)
	finalVitals<-c(timePassed=timePassed,finalVitals)
		#time
		#okayTime<-(timePassed>=runConditions[[1]][1] & timePassed<=runConditions[[1]][2])
	#other vitals
	okayVitals<-sapply(1:4,function(i){
		var<-finalVitals[i]
		varRange<-runConditions[[i]]
		var>=varRange[1] & var<=varRange[2]
		})
	#finalCheck<-all(finalVitals)
	if(any(!okayVitals)){
		print(finalVitals)
		#browser()
		stop(paste0("Accepted run as outside of bounds set for conditions:",
			names(runConditions)[!okayVitals],collapse=", "))
		}
	finalCheck<-all(okayVitals)
	return(finalVitals)
	}		
	
