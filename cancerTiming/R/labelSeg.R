#function to give label as to what side of centromere

# ##use hg19 build for coordinates of centromeres:
# hg19chromosomes<-read.table("~/Documents/Sequencing/CancerTiming/chromosomeRegions.txt",header=TRUE)
# # hg19chromosomes<-read.table("~/chromosomeRegions.txt",header=TRUE)
# names(hg19chromosomes)[c(1:3)]<-c("chr","start","end")
# hg19chromosomes$start[hg19chromosomes$label=="q"]<-hg19chromosomes$start[hg19chromosomes$label=="q"]+1
# hg19chromosomes$end[hg19chromosomes$label=="p"]<-hg19chromosomes$end[hg19chromosomes$label=="p"]-1
# 
# ##My file, the p is *first* not the *smallest*, so need to fix that
# ##But in fact, this is equivalent for how the chromosomes are ordered
# hg19chromosomes$width<-hg19chromosomes$end-hg19chromosomes$start+1
# hg19chromosomes<-do.call("rbind",by(hg19chromosomes,list(hg19chromosomes$chr),function(x){
# 	whp<-which(x$label=="p")
# 	whq<-which(x$label=="q")
# 	if(x$width[whp]>x$width[whq]){ #switch the labels
# 		x$label[whp]<-"q"
# 		x$label[whq]<-"p"
# 	}
# 	return(x)
# 	
# }))
# row.names(hg19chromosomes)<-NULL
# save(hg19chromosomes,file="~/Documents/RfunctionsGenerally/InternalRPackages/cancerTiming/data/hg19chromosomes.rdata")


#divide each side (p and q) -- only need to do this once and could save as part of the package? Or function that divides
#make the small bits near the centromere and add to near segment... perhaps should discard?
#if(getRversion() >= "2.15.1") utils::globalVariables("hg19chromosomes")
.fpath <- function(name){system.file("extdata", name, package="cancerTiming")}
divideGenome<-function(size=10){
	hg19chromosomes<-read.table(.fpath("hg19chromosomes.txt"),stringsAsFactors=FALSE)
#	load(.fpath("hg19chromosomes.rdata"))
#	data("hg19chromosomes",package="cancerTiming", envir = environment())
	size<-size*1e6
	chr<-hg19chromosomes[hg19chromosomes$label!="centromere",]
	chrSeq<-lapply(1:nrow(chr),function(ii){
		if(chr$label[[ii]]=="p"){
			x<-unlist(chr[ii,c("start","end")])
			start<-seq(x[1],x[2],by=size)
			end<-c(tail(start,-1)+1,x[2])
			
		} 
		else{
			end<-(chr[ii,"end"]-seq(0,chr[ii,"width"],by=size))
			start<-c(tail(end,-1)-1,chr[ii,"start"])
			
		}
		return(data.frame(chr=chr[ii,"chr"],segId=paste(chr[ii,"chr"],chr[ii,"label"],1:length(start),sep=""),start,end))
	})
	return(do.call("rbind",chrSeq))
}

# .mySubjectHits<-function(x){#subjectHits=to
# 	if(packageVersion("IRanges")>="2.5.39"){
# 		S4Vectors::subjectHits(x)
# 	}
# 	else{
# 		IRanges::subjectHits(x)
# 	}
# }
# .myQueryHits<-function(x){#queryHits=from
# 	if(packageVersion("IRanges")>="2.5.39"){
# 	  S4Vectors::queryHits(x)
# 	}
# 	else{
# 		IRanges::queryHits(x)
# 	}
# }

labelSeg<-function(chr,start,end,pctOv=0.10){
	if(requireNamespace("GenomicRanges", quietly = TRUE) & requireNamespace("IRanges", quietly = TRUE)){ #require(GenomicRanges)
		hg19chromosomes<-read.table(.fpath("hg19chromosomes.txt"),stringsAsFactors=FALSE)
		#load(.fpath("hg19chromosomes.rdata"))
#		data("hg19chromosomes",package="cancerTiming", envir = environment())
		if(length(grep("chr",chr))==0) chr<-paste("chr",chr,sep="") #in case its just numbers
		hg19chromosomes<-hg19chromosomes[hg19chromosomes$chr%in%chr,]
		hg19chromosomes$chr<-as.character(hg19chromosomes$chr)


		centGr<-GenomicRanges::GRanges(hg19chromosomes$chr,IRanges::IRanges(hg19chromosomes$start,hg19chromosomes$end),label=hg19chromosomes$label)
		segGr<-GenomicRanges::GRanges(chr,IRanges::IRanges(start,end))
		ov<-GenomicRanges::findOverlaps(subject=segGr,query=centGr)
		#GenomicRanges::subsetByOverlaps(query=segGr,subject=centGr)
		shits <- segGr[S4Vectors::subjectHits(ov)] 
		chits <- centGr[S4Vectors::queryHits(ov)] 
		mint <- GenomicRanges::pintersect(shits, chits)
		spercent <- GenomicRanges::width(mint) / GenomicRanges::width(shits)	
		df<-data.frame(index=S4Vectors::queryHits(ov),pct=spercent,val=GenomicRanges::values(centGr)$label[S4Vectors::queryHits(ov)])
		lab<-by(df,list(S4Vectors::subjectHits(ov)),
			function(x){
				wh<-which(x$pct>=pctOv)
				paste(sort(GenomicRanges::values(centGr)$label[x$index[wh]]),collapse="",sep="")
			})
		labInd<-as.numeric(names(lab)) #just to be safe.
		lab<-as.vector(lab)
		crossCent<-intersect(intersect(grep("centromere",lab),grep("q",lab)),grep("p",lab))
		lab[crossCent]<-"pq"
		whNotCross<-setdiff(grep("centromere",lab),which(lab=="centromere")) #take out those entirely in the centromere
		lab[whNotCross]<-gsub("centromere","",lab[whNotCross])
		lab[lab=="centromere"]<-"c"
		outdf<-data.frame(chr=chr,start=start,end=end)
		outdf$label<-NA
		outdf$label[labInd]<-lab
		return(outdf$label)
	}
  else{cat("Sorry, this little helper function requires the package GenomicRanges >= 1.23.23 and IRanges >= 2.5.39, part of bioConductor. Please install this function before using. Now exiting the function.\n")}
  
}
numChromosome<-function(chr){
	chr<-as.character(chr)
	chr[chr=="X"]<-"23"
	chr[chr=="Y"]<-"24"
	chrN<-as.numeric(chr)
	return(chrN)
}