fcauc.fptp <-
function(FP, TP, TDR, FDR.cut=0.2,do.plot=TRUE)
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
	FDR.sel<-min(which(FDR<=FDR.cut))
	sel<-FDR.sel:l
	
#### find AUC
	
	this.x<-1-FP[sel]
	this.x<-this.x[1:(length(this.x)-1)]+diff(this.x)/2
	this.y<- -diff(TP[sel])
	fcauc<-sum(this.y*this.x)
	
#### generate plot
	
	if(do.plot)
	{
		plot(FP, TP, type="l", col="red", lwd=1.5, xlim=c(0,1), ylim=c(0,1), main=paste("FC.AUC:", signif(fcauc,3), "at FDR", FDR.cut),xlab="false positive rate", ylab="true positive rate")
		lines(FP[sel], TP[sel], col="blue", lwd=2.5)
	
		TP.2<-approx(x=TP[sel], y=FP[sel], xout=seq(0, TP[FDR.sel], length.out=round(TP[FDR.sel]/0.01)), method="linear", ties=max)
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
					points(this.FP, this.TP, col="blue", cex=.1)
				}
			}
		}else{
			message("blue area too small to add to the plot.")
		}
	}
	fcauc
}
