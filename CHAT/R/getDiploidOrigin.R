getDiploidOrigin<-function(dd,oo,para){
	minSQ<-99999
	best.y<-NA
	yy<-seq(-1,-0.05,by=0.05)
	nb.list<-c(0,1,2,3)
	nt.list<-c(1,2,3,4,5,6)
	oo$x<-0
	for(tmp.y in yy){
		foo<-function(i){
			x<-dd[i,6]
			y<-dd[i,4]
			nb<-dd[i,9]
			nt<-dd[i,10]
			if(is.na(nb)|is.na(nt))return(NA)
			x0<-oo$x0
			y0<-tmp.y
			if(nb==1&nt==2){
				delta <- y-y0
				if(delta<0)return(1*abs(delta)) else return(delta)
			}
			mD<-getDistToPath(x0,y0,x,y,type=1,nb,nt)
			return(mD$dist)
		}
		if(!para$is.multicore)SQ<-sum(unlist(lapply(1:nrow(dd),foo)),na.rm=TRUE)
        else{
            SQ<-sum(unlist(mclapply(1:nrow(dd),foo)),na.rm=TRUE)
        }
		if(SQ<minSQ){
			minSQ<-SQ
			best.y<-tmp.y
			#plot(dd[,6],dd[,4],xlim=c(0,0.5),ylim=c(-1,1),col=dd[,10],pch=19,xlab='BAF',ylab='LRR')
			#for(nb in nb.list)
			#for(nt in nt.list){
			#	if(nb>nt/2)next
			#	lines(getCoord(seq(0.1,1,by=0.02),oo$x,best.y,1,nb,nt),col=gray(0.1),lwd=1,lty=2)
			#}
		}
		cat(tmp.y,SQ,'\n')
	}
	foo1<-function(i){
		x<-dd[i,6]
		y<-dd[i,4]
		return(abs(x-oo$x)<para$std.BAF&abs(y-best.y)<para$std.LRR)
	}
	vv<-which(unlist(lapply(1:nrow(dd),foo1)))
	vv12<-which(dd[,4]<=best.y&dd[,9]==1&dd[,10]==2)
	vv<-union(vv,vv12)
	return(oo.new=list(x=oo$x0,y=best.y,list=vv))
}