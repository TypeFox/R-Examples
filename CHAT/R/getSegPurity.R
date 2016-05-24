getSegPurity <- function(seg.dat,oo,AGP=1,type=1,para,rm.thr=50,ref.dd=NULL){
	#Main Inference Function
	# if ref.dd is given, Nb and Nt will not be estimated, but use what is in ref.dd
	#LRR correction
	if(para$is.LRRcorrection){
		seg.dat[seg.dat[,4]>=0,4]<-seg.dat[seg.dat[,4]>=0,4]/para$LRR_correction_amp
		seg.dat[seg.dat[,4]<0,4]<-seg.dat[seg.dat[,4]<0,4]/para$LRR_correction_del
	}
	
	#Remove small segments
	#seg.dat=seg.dat[seg.dat[,7]>=rm.thr,]
	Nseg<-nrow(seg.dat)
    
	
	pointer.list<-rep(0,Nseg)
	pointer.list[oo$list]<-1
	new.sam<-cbind(seg.dat,seg_purity=rep(NA,Nseg),nb=rep(NA,Nseg),nt=rep(NA,Nseg))
	new.sam[oo$list,8]<-0
	new.sam[oo$list,9]<-type
	new.sam[oo$list,10]<-2*type
	#reset BAF and LRR threshold for segment assignments
	#para$std.BAF<-0.02
	#para$std.LRR<-0.08
    
	#Infer segment purity
	nt.list<-1:(2*type+4)
	#png(file<-paste(para$pngdir,rownames(seg.dat)[1],'.png',sep=''))
	cat('Analyzing segmental AGP.')
	if(para$is.plot){
		par(mar=c(4,4,3,1))
		plot(seg.dat[,6],seg.dat[,4],cex=0.1,pch=19,xlim=c(0,0.5),ylim=c(-2,2)+oo$y,xlab='BAF',ylab='LRR',lwd=2,cex.lab=1.5,cex.axis=1.5)
		abline(v=c(0,0.5),col=2,lwd=2,lty=2)
		abline(h=0,col=2,lwd=2,lty=2)
		points(seg.dat[oo$list,6],seg.dat[oo$list,4],pch='+',col='yellow')
	}
	foo<- function(k){
		cat('.',sep='')
		x<-seg.dat[k,6]
		y<-seg.dat[k,4]
		if(pointer.list[k]==1)return(c(0,type,2*type))
		seg.p<-NA
		seg.nb<-NA
		seg.nt<-NA
		min.distance<-999
		if(!is.null(ref.dd)){
			nb<-ref.dd[k,9]
			nt<-ref.dd[k,10]
			if(is.na(nb)|is.na(nt))return(c(NA,NA,NA))
			minDist<-getDistToPath(oo$x,oo$y,x,y,strictness=para$strictness,type=type,nb,nt,is.scaling=para$is.scale)
			p0<-minDist$p.min
			#print(paste(p0,ref.dd[k,8]))
			return(c(p0,nb,nt))
		}
		for(nt in nt.list){
			for(nb in 0:nt){
				if(nb>nt/2)next
				if(nb==type&nt==2*type)next
				if(para$is.plot)lines(getCoord(seq(0.1,1,by=0.02),oo$x,oo$y,type,nb,nt),col=gray(0.1),lwd=1,lty=2)
				minDist<-getDistToPath(oo$x,oo$y,x,y,strictness=para$strictness,type=type,nb,nt,is.scaling=para$is.scale)
				p0<-minDist$p.min
				r0<-1
				if(abs(minDist$x.min-x)<=para$std.BAF*r0&abs(minDist$y.min-y)<=para$std.LRR*r0){
					if(is.na(seg.p)){
						seg.p<-p0
						seg.nb<-nb
						seg.nt<-nt
						min.distance<-minDist$dist
					}
					else{
						#Handle Unidentification Problem: nb=nn or nb=nt-nn
						#Using the mixing ratio closest to AGP
						## minimize objective function:	F= nt-2*type + K*abs(seg.p-AGP)
						F0<- seg.nt-2*type + para$K*abs(seg.p-AGP)
						F1<- nt-2*type + para$K*abs(p0-AGP)
						if(F1<F0){
							seg.p<-p0
							seg.nb<-nb
							seg.nt<-nt
							min.distance<-minDist$dist
						}
					}
				}
			}
		}
		if(para$is.plot)points(seg.dat[k,6],seg.dat[k,4],pch=19,cex=seg.p/AGP,col=seg.nt)
		return(c(seg.p,seg.nb,seg.nt))
	}
	if(para$is.multicore){
        tmp<-mclapply(1:Nseg,foo)
    } else tmp<-lapply(1:Nseg,foo)
	tmp<-t(matrix(unlist(tmp),nrow=3))
	if(para$is.plot){
		points(seg.dat[,6],seg.dat[,4],pch=19,cex=tmp[,1]/AGP,col=tmp[,3])
		points(oo$x,oo$y,pch='*',cex=2,col='purple')
	}
	new.sam[,8:10]<-tmp
	if(para$is.plot)title(paste(rownames(seg.dat)[1],'Ploidy:',type*2))
	#dev.off()
	return(new.sam)
}