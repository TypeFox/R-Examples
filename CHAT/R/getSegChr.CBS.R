getSegChr.CBS <- function(bb.chr=NULL,ll.chr=NULL,sam.col=5,control=TRUE,thr.hets=0.1,data.type='copy',controlOne=0,nt=FALSE){
	##Process segmentation by chromosome.
	chr<-bb.chr[1,2]
	cat('processing ',chr,'\n')
	if(nrow(bb.chr)!=nrow(ll.chr))stop("Incompatible markers for BAF and LogR!")
	nsam<-ncol(bb.chr)
	bb.chr<-bb.chr[order(as.numeric(bb.chr[,3])),]
	ll.chr<-ll.chr[order(as.numeric(ll.chr[,3])),]
	dat.chr<-c()
	i<-sam.col
    if(control)delta=2 else delta=1
    if(controlOne>0)delta=1
	sampleID<-colnames(bb.chr)[i:nsam]
    foo<-function(i){
		id<-colnames(bb.chr)[i]
		paired.id<-colnames(bb.chr)[i+delta-1]
        if(i+delta-1>nsam)stop("Check your data. Are all samples paired?")
        if(controlOne>0)paired.id<-colnames(bb.chr)[sam.col+controlOne-1]
		tag<-getHets(baf=bb.chr[,id],paired.baf=bb.chr[,paired.id],thr=thr.hets)
		pp.bb<-bb.chr[tag,3]
		pp.ll<-ll.chr[,3]
		bb<-bb.chr[tag,id]
		ll<-ll.chr[,id]
		ll.paired<-ll.chr[,paired.id]
		CNA.bb<-CNA(abs(bb-0.5),bb.chr[tag,2],pp.bb,sampleid=id)
		s.CNA.bb<-smooth.CNA(CNA.bb)
		seg.b<-segment(s.CNA.bb,min.width=5)
		CNA.ll<-CNA(ll,ll.chr[,2],pp.ll,sampleid=id)
		s.CNA.ll<-smooth.CNA(CNA.ll)
		seg.l<-segment(s.CNA.ll,min.width=5)
		seg<-MergeBreakPointsByChr(chr=chr,id=id,seg.b$output,seg.l$output,bb.chr[,3])
		ss<-seg[,3:4]
		nsegs<-nrow(seg)
		if(is.null(nsegs)){
            i<-i+delta
            next
        }  ## hromosome contains no segments (potential error)
        tmp.dd<-c()
		for(j in 1:nsegs){
			tmp.bb<-which(pp.bb>=as.numeric(seg[j,3])&pp.bb<=as.numeric(seg[j,4]))
			tmp.ll<-which(pp.ll>=as.numeric(seg[j,3])&pp.ll<=as.numeric(seg[j,4]))
			BAF.mean<-getBAFmean(bb[tmp.bb])
			if(length(tmp.ll)==0)logR.mean<-0 else{
				if(data.type=='copy')logR.mean<-log2(median(ll[tmp.ll]/ll.paired[tmp.ll]^0,na.rm=TRUE))-1
				if(data.type=='log')logR.mean=median(ll[tmp.ll],na.rm=TRUE)
            }
			BAF.num<-length(tmp.bb)
			logR.num<-length(tmp.ll)
			sam.seg<-cbind(id,chr,seg[j,3],seg[j,4],logR.mean,logR.num,BAF.mean,BAF.num)
            tmp.dd <- rbind(tmp.dd,sam.seg)
			cat('.')
		}
        #i<-i+delta
        #if(i>nsam)break
		cat('\ndone!\n')
        return(tmp.dd)
	}
    if(nt&controlOne==0&control){
        tmp=mclapply(seq(i,nsam,by=2),foo)
        dat.chr=c()
        for(kk in tmp)dat.chr=rbind(dat.chr,kk)
    }else{
        dat.chr=c()
        while(i<=nsam){
            sam.seg=foo(i)
            dat.chr=rbind(dat.chr,sam.seg)
            i=i+delta
            if(i>nsam)break
        }
    }
	return(na.omit(dat.chr))
}
