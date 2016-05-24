#######analyze each of the data sets
eventTimingOverList<-function(dfList,normCont,eventArgs=NULL){
	###used for timing straightforward 
	#assume dfList is list of data per sample; each sample must have data frame consisting of all mutations with the following columns:
		#type one of c("Other","CNLOH","SingleGain","Diploid","DoubleGain")
		#segId saying for each mutation, what segment it is in. Can not contain a "."
		#nMutAllele == number of reads with mutation
		#nReads == total number of reads
		#mutationId -- a unique id for each mutation
	#eventArgs=list of arguments passed to eventTiming via 'do.call'. Should NOT contain the arguments 'x','m','history','totalCopy','type','mutationId' or 'normCont'

	if(any(sapply(dfList,function(x){any(!c("segId","type","nMutAllele","nReads","mutationId")%in%names(x))}))) stop("dfList elements are missing required column names")
	if(any(sapply(dfList,function(x){length(grep("[.]",x$segId))>0}))) stop("segId cannot contain a period")
	if(is.null(names(dfList))) names(dfList)<-paste("Sample",1:length(dfList),sep="")
	if(length(normCont)!=length(dfList)) stop("dfList and normCont must be of the same length (1 per sample)")

	whCall<-lapply(dfList,function(x){tapply(x$segId,factor(x$type,levels=c("Other","CNLOH","SingleGain","Diploid","DoubleGain")),unique,simplify = FALSE)})
	singleSampleFunction<-function(x,nc,dat){
		if(length(x[["CNLOH"]])>0){
			ACNLOH<-makeEventHistory(totalCopy=2,type="LOH")[[1]]
			eventCNLOH<-lapply(x[["CNLOH"]],function(segId){
				subdat<-dat[which(dat$segId==segId),]
				#print(segId)
		#		print(dim(subdat))
				out<-do.call(eventTiming,c(list(x=subdat$nMutAllele, m=subdat$nReads, 
					history=ACNLOH,totalCopy=2,type="CNLOH",mutationId=subdat$mutationId,
					normCont=nc),eventArgs))
			})
			names(eventCNLOH)<-x[["CNLOH"]]
		}
		else eventCNLOH<-NULL
		if(length(x[["SingleGain"]])>0){
			AGain<-makeEventHistory(totalCopy=3,type="gain")[[1]]
			eventGain<-lapply(x[["SingleGain"]],function(segId){
				subdat<-dat[which(dat$segId==segId),]
				#print(segId)
				# print(dim(subdat))
				do.call(eventTiming,c(list(x=subdat$nMutAllele, m=subdat$nReads, 
					history=AGain,totalCopy=3,type="gain",mutationId=subdat$mutationId,
					normCont=nc),eventArgs))
			})
			names(eventGain)<-x[["SingleGain"]]			
		}
		else eventGain<-NULL
		if(length(x[["DoubleGain"]])>0){
			ADGain<-makeEventHistory(totalCopy=4,type="gain")[[1]]
			eventDGain<-lapply(x[["DoubleGain"]],function(segId){
			subdat<-dat[which(dat$segId==segId),]
			#print(segId)
			# print(dim(subdat))
			do.call(eventTiming,c(list(x=subdat$nMutAllele, m=subdat$nReads, 
				history=ADGain,totalCopy=4,type="gain",mutationId=subdat$mutationId,
				normCont=nc),eventArgs))
				})
			names(eventDGain)<-x[["DoubleGain"]]
		}
		else eventDGain<-NULL
		return(list(SingleGain=eventGain,CNLOH=eventCNLOH,DoubleGain=eventDGain))
	}
	mapply(whCall,normCont,dfList,FUN=singleSampleFunction,SIMPLIFY=FALSE)
}

getPi0Summary<-function(eventList,CI=TRUE){ 
	xmle<-lapply(eventList,.getPiSingle,CI=CI)
	xmle<-data.frame(Sample=rep(names(xmle),times=sapply(xmle,nrow)),do.call(rbind,xmle))
	rankWInSample<-do.call("rbind",tapply(1:nrow(xmle),xmle$Sample, function(ii){
		p<-xmle$pi0[ii]
		r<-rep(NA,length(p))
		r[!is.na(p)]<-rank(p[!is.na(p)])
		return(data.frame(rank=r,segId=xmle$segId[ii],Sample=xmle$Sample[ii]))
	}))
	xmle$rankInSample<-rankWInSample$rank[match(paste(xmle$Sample,xmle$segId),paste(rankWInSample$Sample,rankWInSample$segId))]
	return(xmle)
}

.getPiSingle<-function(estList,CI){
	pi0<-sapply(unlist(estList,recursive=FALSE),function(x){x$pi["Stage0"]})
	if(length(pi0)>0){
		vals<-strsplit(names(unlist(estList,recursive=FALSE)),"[.]")	
		nam<-sapply(vals,.subset2,2)
		type<-sapply(vals,.subset2,1)
		N<-sapply(unlist(estList,recursive=FALSE),function(x){x$summaryTable[2,]})
		if(CI){
			ui<-sapply(unlist(estList,recursive=FALSE),function(x){if("piCI" %in% names(x)) x$piCI["Stage0",2] else NA})
			li<-sapply(unlist(estList,recursive=FALSE),function(x){if("piCI" %in% names(x)) x$piCI["Stage0",1] else NA})
			out<-data.frame(pi0=pi0,lCI=li,uCI=ui,N=N,type=type,segId=nam)			
		}
		else{
			out<-data.frame(pi0=pi0,N=N,type=type,segId=nam)			
		}
		row.names(out)<-NULL
		out<-out[order(out$segId),]
		return(out)
	}
	else{
		if(CI) return(data.frame(pi0=NA,lCI=NA,uCI=NA,N=NA,type=NA,segId=NA))
		else return(data.frame(pi0=NA,N=NA,type=NA,segId=NA))
	}
}