getSampleAGP <- function(sam.dat,oo,para){
	if(oo$x>=para$thr.triploidy & oo$x<=0.167)para$is.tri<-TRUE
	else{para$is.tri<-FALSE}
	Ns<-dim(sam.dat)[1]
	xx<-1:Ns
	kk<-sam.dat[!xx%in%oo$list,5]
	QC.GC<-ifelse(length(kk)>0,sum(kk)/sum(sam.dat[,5]),0)
	p1.list<-seq(0.05,0.95,by=0.1)
	SD.min<-99999
	cali.best<-NULL
	id<-rownames(sam.dat)[1]
	typeList<-ifelse(para$is.tri,1.5,list(1:4))
	modelList<-ifelse(para$is.tri,1,list(2))
	SD1<-999999
	SD2<-999999
	for(p in p1.list){
		for(tp in unlist(typeList)){
			for(model in unlist(modelList)){
				para$model<-model
				cali<-getSumDist(oo,p,sam.dat,type=tp,para=para)
				if(cali$SD<SD.min){
					SD.min<-cali$SD
					QC.CL<-sum(sam.dat[cali$clist,5])/(sum(sam.dat[,5])-sum(sam.dat[oo$list,5]))
					QC.CP<-sum(sam.dat[cali$plist,5])/(sum(sam.dat[cali$clist,5]))
					cali.best<-cali
					cali.best$type<-tp
					cali.best$sam.p<-p
					cali.best$model<-model
					cali.best$percent.change<-QC.GC
					cali.best$percent.on.track<-QC.CL
					cali.best$percent.on.point<-QC.CP
					#png(file=paste('D://bo/tmp/',p,tp,'.png'))
					if(para$is.plot)plotBAFLRR(sam.dat,cali.best,oo,para=para)
					#dev.off()
				}
				if(model==1 & cali$SD<SD1)SD1<-cali$SD
				if(model==2 & cali$SD<SD2)SD2<-cali$SD
			}
		}
	}
	p2.list<-seq(max(cali.best$sam.p-0.1,0.01),min(cali.best$sam.p+0.1,0.99),by=0.01)
	for(p in p2.list){
		for(tp in unlist(typeList)){
			for(model in unlist(modelList)){
				para$model<-model
				cali<-getSumDist(oo,p,sam.dat,type=tp,para=para)
				if(cali$SD<SD.min){
					SD.min<-cali$SD
					QC.CL<-sum(sam.dat[cali$clist,5])/(sum(sam.dat[,5])-sum(sam.dat[oo$list,5]))
					QC.CP<-sum(sam.dat[cali$plist,5])/(sum(sam.dat[cali$clist,5]))
					cali.best<-cali
					cali.best$type<-tp
					cali.best$sam.p<-p
					cali.best$model<-model
					cali.best$percent.change<-QC.GC
					cali.best$percent.on.track<-QC.CL
					cali.best$percent.on.point<-QC.CP
					#png(file=paste('D://bo/tmp/',p,tp,'.png'))
					if(para$is.plot)plotBAFLRR(sam.dat,cali.best,oo,para=para)
					#dev.off()
				}
				if(model==1 & cali$SD<SD1)SD1<-cali$SD
				if(model==2 & cali$SD<SD2)SD2<-cali$SD
			}
		}
	}
	para$model<-cali.best$model
	cali.best$likelihood<-c(SD1,SD2)
	if(para$is.png){
		png(filename=paste(para$pngdir,'/',id,'.png',sep=''))
		plotBAFLRR(sam.dat,cali.best,oo,para=para)
		dev.off()
	}
	if(para$is.perm){
		plist<-getPermutation(sam.dat,oo,type=cali.best$type,para=para)
		if(is.na(plist)){conf.int<-c(NA,NA,NA)}
		else{conf.int<-quantile(plist,c(0.025,0.5,0.975))}
		cali.best$conf.int<-conf.int
		cali.best$std<-sqrt(var(plist,na.rm=TRUE))
		print(cali.best$std)
	}
	return(cali.best)
}
