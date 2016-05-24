getAGP <- function(para){
	load(para$datafile)
	dd.dat<-dd.dat[dd.dat[,7]>=para$BAFfilter,]
	if(para$is.normalize)dd.dat<-NormalizeLRR(dd.dat,para=para)
	sID<-unique(rownames(dd.dat))
	if(para$save.origin)OO<-c()
	if(para$is.perm)write('sampleid\tpurity\tstd\tPC\tPoP\tEuploidy\tpl.tumor\tpl.overall\t%amp\t%del\t%del.loh\t%cn.loh\t0.025%\tmedian\t0.975%',file=para$savefile)
	else	write('sampleid\tpurity\tPC\tPoP\tEuploidy\tpl.tumor\tpl.overall\t%amp\t%del\t%del.loh\t%cn.loh',file=para$savefile)
	for(id in sID){
		cat('processing ',id)
		sam.dat<-dd.dat[rownames(dd.dat)==id,]
		sam.dat<-na.omit(sam.dat)
		sam.dat<-sam.dat[is.finite(sam.dat[,4]),]
		print(nrow(sam.dat))
		if(length(sam.dat)==0){next}
		#sam.dat=sam.dat[order(sam.dat[,6]),]
		if(!is.null(para$exclude.chr)){
			sam.dat<-sam.dat[which(!sam.dat[,1]%in%para$exclude.chr),]
		}
		if(nrow(sam.dat)>=1000){
			tmp.thr<-sort(sam.dat[,7],decreasing=TRUE)[1000]
			sam.dat<-sam.dat[sam.dat[,7]>=tmp.thr,]
		}
		oo<-getOrigin(sam.dat,para=para)
		if(para$save.origin)OO<-c(OO,list(oo))
		cali<-getSampleAGP(sam.dat,oo,para=para)
		cali<-getLOH(cali,sam.dat,para=para)
		AmpDel<-getAmpDel(oo,sam.dat)
		Ploidy<-getPloidy(cali$sam.p,cali$type,sam.dat,oo=oo)
		cat(' done\n')
		type<-cali$type
		nn<-2*type
		del.loh<-cali$del.loh
		cn.loh<-cali$cn.loh
		dl<-0
		cl<-0
		if(length(del.loh)!=0){dl<-sum(sam.dat[del.loh,5])/sum(sam.dat[,5])}
		if(length(cn.loh)!=0){cl<-sum(sam.dat[cn.loh,5])/sum(sam.dat[,5])}
		ci<-cali$conf.int
		model<-cali$model
		if(para$is.perm){
			write(paste(id,cali$sam.p,signif(cali$std,digits=3),signif(cali$percent.change,digits=3),signif(cali$percent.on.point*cali$percent.on.track,digits=3),nn,signif(Ploidy[2],digits=3),signif(Ploidy[1],digits=3),signif(AmpDel[1],digits=3),signif(AmpDel[2],digits=3),signif(dl,digits=3),signif(cl,digits=3),ci[1],ci[2],ci[3],sep='\t'),file=para$savefile,append=T)
		}
		else{
			write(paste(id,cali$sam.p,signif(cali$percent.change,digits=3),signif(cali$percent.on.point*cali$percent.on.track,digits=3),nn,signif(Ploidy[2],digits=3),signif(Ploidy[1],digits=3),signif(AmpDel[1],digits=3),signif(AmpDel[2],digits=3),signif(dl,digits=3),signif(cl,digits=3),sep='\t'),file=para$savefile,append=T)
		}
	}
	if(para$save.origin){
		names(OO)<-sID
		save(OO,file='Origin.Rdata')
	}
}	## Main function: AGP estimation