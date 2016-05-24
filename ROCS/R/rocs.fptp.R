rocs.fptp <-
function(FP, TP, TDR, FDR.cut=0.2)
{
	
	o<-order(-FP, -TP)
	FP<-FP[o]
	TP<-TP[o]
	TDR<-TDR[o]
	
	FP<-c(1, FP, 0)
	TP<-c(1, TP, 0)
	TDR<-c(min(TDR), TDR, 1)
	FDR<-1-TDR
	l<-length(FP)
	
	abs.d.TP<- -diff(TP)
	mid.FP<- FP[1:(l-1)]+diff(FP)/2	
	mid.TDR<-TDR[1:(l-1)]+diff(TDR)/2
	
	vus<-sum(abs.d.TP*(1-mid.FP)*mid.TDR)
	
	clear3d()
	lines3d(FP, TP, TDR, col="red", lwd=2)
	lines3d(FP, TP, rep(0,l),col="red",lwd=1.5)
	colorlut <- topo.colors(100)
	col<-colorlut[round(mid.TDR/max(TDR)*99)+1]
	
	r<-matrix(0, nrow=6*(l-1), ncol=3)
	r.pointer<-1
	
	for(i in 1:(l-1)) 
	{
		r[r.pointer:(r.pointer+2),]<-cbind(c(FP[i], FP[i+1],1), c(TP[i], TP[i+1], TP[i]), c(TDR[i], TDR[i+1],TDR[i]))
		r[(r.pointer+3):(r.pointer+5),]<-cbind(c(1,1,FP[i+1]), c(TP[i], TP[i+1], TP[i+1]), c(TDR[i], TDR[i+1],TDR[i+1]))
		
		r.pointer<-r.pointer+6
	}
	triangles3d(r,col=rep(col, each=6))
	
	
	decorate3d(aspect="iso",xlab="FPR",ylab="TPR",zlab="TDR",xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), box=TRUE, axes=TRUE)
	aspect3d(1,1,1)
	
	FDR<- -cummax(-FDR)
	l<-length(TP)
	FDR.sel<-min(which(FDR<=FDR.cut))
	sel<-FDR.sel:l
	
	
	this.x<-1-FP[sel]
	this.x<-this.x[1:(length(this.x)-1)]+diff(this.x)
	this.y<- -diff(TP[sel])
	fcauc<-sum(this.y*this.x)
	
	lines3d(x=c(FP[FDR.sel],1), y=c(TP[FDR.sel],TP[FDR.sel]), z=c(0,0))
	
	TP.2<-approx(x=TP[sel], y=FP[sel], xout=seq(0, TP[FDR.sel], length.out=round(TP[FDR.sel]/0.01)), method="constant", ties=max)
	if(length(TP.2$x)>0)
	{
		for(i in 1:length(TP.2$x))
		{
			all.FP<-seq(0,1,by=0.01)
			if(!is.na(TP.2$y[i]))
			{
				this.FP<-all.FP[all.FP>=TP.2$y[i]]
				if(length(this.FP)>0)
				this.TP<-rep(TP.2$x[i], length(this.FP))
				points3d(this.FP, this.TP, rep(0, length(this.TP)), col="grey", cex=.1)
			}
		}
	}else{
		message("FDR-restricted AUC area too small to add to the plot.")
	}
	
	
	message(paste("At the FDR cutoff of", FDR.cut, ", the restricted AUC (shade) is:", fcauc))
	
	vus	
}
