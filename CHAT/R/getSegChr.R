getSegChr <- function(bb.chr=NULL,ll.chr=NULL,sam.col=5,control=TRUE, thr.hets=0.1, bin=1000,data.type='copy',controlOne=0){
	chr<-bb.chr[1,2]
	cat('processing ',chr,'\n')
	#desired input for bb.chr and ll.chr:
	##row_names, col1, col2, col3,... coln
	## SNP id, chr, pos, sample1,... samplen
	## if is.paired is TRUE, for TCGA samples only, i.e. when is.paired=T, is.tcga is T
	## if is.uniform is TRUE, then all samples share same segments, and the first column of segment file will not be used
	
	if(nrow(bb.chr)!=nrow(ll.chr))stop("Incompatible markers for BAF and LogR!")
	nsam<-ncol(bb.chr)
	bb.chr<-bb.chr[order(as.numeric(bb.chr[,3])),]
	ll.chr<-ll.chr[order(as.numeric(ll.chr[,3])),]
	#bb.chr=bb.chr[order(bb.chr[,2]),]
	dat.chr<-c()
	i<-sam.col
	sampleID<-colnames(bb.chr)[i:nsam]
    if(control)delta<-2 else delta<-1
    if(controlOne>0)delta<-1
	while(i<=nsam){
		id<-colnames(bb.chr)[i]
		paired.id<-colnames(bb.chr)[i+delta-1]
        if(i+delta-1>nsam)stop("Check your data. Are all samples paired?")
        if(controlOne>0)paired.id<-colnames(bb.chr)[sam.col+controlOne-1]
		tag<-getHets(baf=bb.chr[,id],paired.baf=bb.chr[,paired.id],thr=thr.hets)
		pp.bb<-bb.chr[,3]
		pp.ll<-ll.chr[,3]
		bb<-bb.chr[,id]
		ll<-ll.chr[,id]
        ll.paired<-ll.chr[,paired.id]
        if(id==paired.id)A<-0 else A<-0
		Nt<-length(tag)
		kk<-1
		while(kk<Nt){
			ed<-kk+bin-1
			if(ed>Nt)ed<-Nt
			vv<-tag[kk:ed]
			BAF.mean<-getBAFmean(bb[vv])
			if(data.type=='copy')logR.mean<-log2(median(ll[vv]/ll.paired[vv]^A, na.rm=TRUE))-1
            if(data.type=='log')logR.mean<-median(ll[vv]-ll.paired[vv]*A,na.rm=TRUE)
			sam.seg<-cbind(id,chr,pp.bb[vv[1]],pp.bb[vv[length(vv)]],logR.mean,bin,BAF.mean,bin)
			dat.chr<-rbind(dat.chr,sam.seg)
			kk<-ed+1
			cat('.')
		}
        i<-i+delta
        if(i>nsam)break
		cat('\ndone!\n')
	}
	return(na.omit(dat.chr))
}