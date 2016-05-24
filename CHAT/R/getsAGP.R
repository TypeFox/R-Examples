getsAGP <- function(para){
	## Main Function 2
	new.dd<-c()
	cat('loading files.')
	load(para$inputdata)
    if(para$is.normalize)dd.dat=NormalizeLRR(dd.dat,para=para)
	cat('done\n')
	#load(para$originFile)
	pG<-as.matrix(read.table(para$purityfile,header=T,stringsAsFactors=F,sep='\t'))
	rownames(pG)<-pG[,1]
	sID<-unique(rownames(dd.dat))
	for(i in 1:length(sID)){
		id<-sID[i]
		pl<-as.numeric(pG[id,5])/2
		cat('processing ',id,pl,'\t')
		seg.dat<-unique(dd.dat[rownames(dd.dat)==id,])
		seg.dat<-na.omit(seg.dat)
		seg.dat<-seg.dat[is.finite(seg.dat[,4]),]
		if(nrow(seg.dat)>=1000){
			tmp.thr<-sort(seg.dat[,7],decreasing=T)[1000]
			tmp.dat<-seg.dat[seg.dat[,7]>=tmp.thr,]
			posID<-paste(tmp.dat[,1],tmp.dat[,2],tmp.dat[,3],sep='-')
			rownames(seg.dat)<-paste(seg.dat[,1],seg.dat[,2],seg.dat[,3],sep='-')
			oo<-getOrigin(tmp.dat,para=para)
			oo$list<-match(posID[oo$list],rownames(seg.dat))
			rownames(seg.dat)<-rep(id,nrow(seg.dat))
			dd<-getSegPurity(seg.dat,oo,AGP=as.numeric(pG[id,2]),para=para,type=pl)
		}else{
			oo<-getOrigin(seg.dat,para=para)
			dd<-getSegPurity(seg.dat,oo,AGP=as.numeric(pG[id,2]),para=para,type=pl)
		}
		if(pl>=2){
            ## recalculate tetraploid tumors
            oo.new<-getDiploidOrigin(dd,oo,para=para)
			#para$is.plot=T
			#png(file=paste('tetraBRCA-',id,'.png',sep=''))
			vv0<-which(dd[,9]==1&dd[,10]==2)
			dd[vv0,9]<-pl
			dd[vv0,10]<-2*pl
			dd<-getSegPurity(dd[,1:7],oo.new,AGP=as.numeric(pG[id,2]),para=para,type=1,ref.dd=dd)
			#dev.off()
			#para$is.plot=F
		}
		if(para$sensitivity){
			cat('Performing sensitivity analysis.\n')
			SS<-NULL
			oo.x<-runif(para$ss.N,max(oo$x-para$ss.BAF,0),min(oo$x+para$ss.BAF,0.5))
			oo.y<-runif(para$ss.N,oo$y-para$ss.LRR,oo$y+para$ss.LRR)
			SS$results<-vector(para$ss.N,mode='list')
			for(i in 1:para$ss.N){
				cat('.')
				tmp.oo<-oo
				tmp.oo$x<-oo.x[i]
				tmp.oo$y<-oo.y[i]
				tmp.dd<-getSegPurity(seg.dat,tmp.oo,AGP=as.numeric(pG[id,2],para=para,type=as.numeric(pG[id,5]/2)))
				SS$results[[i]]<-tmp.dd[,c(1:3,8:10)]
			}
			SS$oo<-oo
			SS$x<-oo.x
			SS$y<-oo.y
			SS$seg.p<-dd[,c(1:3,8:10)]
			save(SS,file=paste(para$ss.dir,id,'-ss.Rdata',sep=''))
		}
		new.dd<-rbind(new.dd,dd)
		cat('done\n')
	}
	save(new.dd,para,file=para$savedata)
    
}	## Main function: sAGP estimation