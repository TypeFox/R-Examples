

mut2Seg<-function(mutData,segData,verbose=TRUE){
  if(requireNamespace("GenomicRanges", quietly = TRUE) && requireNamespace("IRanges", quietly = TRUE)){ #require(GenomicRanges)
    
		if(any(!c('chr','end','start') %in% colnames(segData))) stop("segData must have columns with names 'chr','start',and 'end'")
		else{
			segData<-data.frame(segData)
			if(!'strand' %in% colnames(segData)) segData$strand<-"*"
			if(!'width' %in% colnames(segData)) segData$width<-segData$end-segData$start+1
			wh<-sapply(c('chr','start','end','strand','width'),grep,colnames(segData))
			if(length(wh)<ncol(segData)) segData<-cbind(segData[,wh],segData[,-wh])
			else segData<-segData[,wh]
		}
		if(!all(c("chr","position")%in%colnames(mutData))) stop("column names of mutData must include 'chr' and 'position'")
	
		#check same chromosome names:
		segChr<-unique(segData$chr)
		mutChr<-unique(mutData$chr)
		newChr<-union(segChr,mutChr)
		if(!all(newChr%in%intersect(segChr,mutChr))){
			#warning("Not all of the chromosomes shared in mutData and segData, will allow the union of them")
			if(any(!newChr%in%segChr)){
				if(verbose){
					cat("Chromosomes not in segData:\n")
					print(newChr[!newChr%in%segChr])
				}
			}
			if(any(!newChr%in%mutChr)){
				if(verbose){
					cat("Chromosomes not in mutData:\n")
					print(newChr[!newChr%in%mutChr])
				}
			}
		}	
		segGr<-GenomicRanges::GRanges(seqnames =factor(segData$chr,levels=newChr),
		          ranges = IRanges::IRanges(segData$start, end=segData$end),strand="*",
		          segData[,6:ncol(segData)])
		mutGr<-GenomicRanges::GRanges(seqnames=factor(mutData$chr, levels=newChr), ranges=IRanges::IRanges(mutData$position,end=mutData$position))
		ov<-GenomicRanges::findOverlaps(segGr,mutGr)
		qHits<-S4Vectors::queryHits(ov) #queryHits=from
		sHits<-S4Vectors::subjectHits(ov) #subjectHits=to
		combDf<-data.frame(mutData[sHits,],
			seg_start=start(segGr)[qHits],seg_end=end(segGr)[qHits],
				as.data.frame(GenomicRanges::values(segGr)[qHits,]))

		#check that none match more than 1
		if(length(unique(sHits))!=length(sHits)){
			isDup<-duplicated(sHits)
			whDup<-which(sHits%in%sHits[which(isDup)])
			ndups<-table(sHits[whDup])
			if(verbose){
				cat("Overlapping Segments with Mutations matching:\n")
				print(segData[unique(qHits[whDup]),])
			}
	#		stop(length(ndups)," mutations matched more than 1 segment")
		
		}
		#add those with no hits!
		nmissing<-nrow(mutData)-length(sHits)
		if(nmissing>0){
			dummyData<-data.frame(seg_start=NA,seg_end=NA,matrix(NA,ncol=ncol(GenomicRanges::values(segGr)),nrow=nmissing))
			names(dummyData)<-c("seg_start","seg_end",colnames(GenomicRanges::values(segGr)))
			if(length(sHits)>0) missDat<-cbind(mutData[-sHits(ov),],dummyData)
			else{
				missDat<-cbind(mutData,dummyData) #means there were no mutations in these segments!
				if(verbose) cat("there was no overlap between the mutations and the segments\n")
			}
			colnames(missDat)<-make.names(colnames(missDat))
			colnames(combDf)<-make.names(colnames(combDf))
			combDf<-rbind(combDf,missDat)
		}
		combDf<-combDf[order(numChromosome(combDf$chr),combDf$position),]
		return(combDf)
	}
  else{cat("Sorry, this little helper function requires the package GenomicRanges and IRanges, part of bioConductor. Please install this function before using. Now exiting the function.\n")}
  
}