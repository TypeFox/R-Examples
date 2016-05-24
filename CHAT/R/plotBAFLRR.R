plotBAFLRR <-
function(sam.dat,cali.best,oo,origin='on',para=para){
	gg<-getGrid(oo$x0,oo$y0,cali.best$sam.p,Nc=8,type=cali.best$type,para=para)
	if(para$is.bubble){
		plot(sam.dat[,6],sam.dat[,4],cex.axis=1.5,cex.lab=1.5,lwd=3,cex=sam.dat[,5]/1000,main=paste('AGP: ',signif(cali.best$sam.p,digits=3),'   Raw Ploidy: ',cali.best$type,'\nPC: ', signif(cali.best$percent.change,digits=3),'  PoP: ',signif(cali.best$percent.on.point*cali.best$percent.on.track,digits=3),sep=''),xlab='BAF',ylab='LRR',xlim=c(0,0.5),ylim=c(-2,2))
	}
	else{
		plot(sam.dat[,6],sam.dat[,4],main=paste('purity: ',signif(cali.best$sam.p,digits=3),'  ploidy:',cali.best$type,sep=''),xlab='BAF',ylab='LRR',xlim=c(0,0.5),ylim=c(-2,2),cex=0.1)
	}
	abline(0,0,col='brown',lwd=2)
	abline(v=0.5,col='brown',lwd=2)
	abline(v=0,col='brown',lwd=2)
	for(i in 1:para$num.tracks){
		#abline(cali.best$b[i],cali.best$k[i],col='blue',lty=2,lwd=2)
	}
	if(origin=='on'){
		points(sam.dat[oo$list,6],sam.dat[oo$list,4],pch=3,cex=0.6,col='yellow',lwd=2)
	}
	p<-cali.best$sam.p
	type<-cali.best$type
	nn<-type
	if(para$model==1){
		n0<-log2(1+(nn-1)*p)
		nn<-1
	}
	else{
		n0<-log2(nn)
	}
	if(para$is.tri){
		n0<-log2(2+p)-1
		b0<-abs(1/(2+p)-0.5)
		nn<-1
	}
	points(gg,pch=8,col='red',lwd=2)
	size.line<-100
	p<-seq(0,1,length.out=size.line)
	for(nt in 1:8){
		ll<-log2(2*nn*(1-p)+nt*p)-1-n0+oo$y0
		for(nb in 0:nt){
			if(nb>nt/2)next
			bb<-abs((nb*p+nn*(1-p))/(2*nn*(1-p)+nt*p)-0.5)+oo$x0
			if(para$is.tri)bb<-abs((p*nb+(1-p))/(p*nt+(1-p)*2*nn)-0.5)+oo$x0-b0
			lines(bb,ll,col='red',lwd=2,lty=2)
		}
	}
}
