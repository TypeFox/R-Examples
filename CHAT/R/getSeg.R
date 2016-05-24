getSeg <- function(BAFfilePath,logRfilePath,sam.col=5,control=TRUE,thr.hets=0.1,bin=500,savefile,data.type='copy', cbs=FALSE, controlOne=0,nt=FALSE){
	
	#Must provide the full path to BAF and logR files
	#Each folder must contain 22-24 files, each for one chromosome. 23 for chrX, 24 for chrY
    
	dd<-c()
	BAFfiles<-dir(BAFfilePath,full.names=T)
	LRRfiles<-dir(logRfilePath,full.names=T)
    Nc<-length(BAFfiles)
    chrlist<-1:Nc
	for(chr in chrlist){
		bb.chr<-read.table(BAFfiles[chr],header=T,stringsAsFactors=F)
		ll.chr<-read.table(LRRfiles[chr],header=T,stringsAsFactors=F)
        if(!cbs)cc<-getSegChr(bb.chr,ll.chr,bin=bin,data.type=data.type,control=control,thr.hets=thr.hets,sam.col=sam.col,controlOne=controlOne) else cc<-getSegChr.CBS(bb.chr,ll.chr,data.type=data.type,control=control,thr.hets=thr.hets,sam.col=sam.col,controlOne=controlOne,nt=nt)
		dd<-rbind(dd,cc)
	}
	dd.dat<-as.matrix(dd[,2:8])
	dd.dat<-matrix(as.numeric(dd.dat),ncol=7)
	colnames(dd.dat)<-colnames(dd)[2:8]
	rownames(dd.dat)<-dd[,1]
	dd.dat<-dd.dat[order(rownames(dd.dat)),]
	save(dd.dat,file=savefile)
}	## Main function: Binned segmentation
