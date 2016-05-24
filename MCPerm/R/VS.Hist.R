VS.Hist <-
function(Trad_data,MC_data,Trad_col="grey",MC_col="black",title=NULL,xlab=NULL){
    ## data=LnOR,VARLnOR,Qp,I2,merged_LnOR,merged_VARLnOR,merged_LnOR_p
	 
	 xlim_value=range(c(Trad_data,MC_data))
	 Trad_ylim=max(hist(Trad_data,plot=FALSE)$density)
	 MC_ylim=max(hist(MC_data,plot=FALSE)$density)
	 ylim_value=c(0,max(Trad_ylim,MC_ylim))
		 
    hist(Trad_data,freq=FALSE,border=Trad_col,
	       main=title,xlim=xlim_value,ylim=ylim_value,col=Trad_col,xlab=xlab)
	 
	 par(new=TRUE)
	 hist(MC_data,freq=FALSE,border=MC_col,
	       main="",xlim=xlim_value,ylim=ylim_value,xlab="")
			 
	 legend("top",legend=c("TradPerm","MCPerm"), fill=c(Trad_col,"white"),
		 border=c("white",MC_col),cex=0.7,xpd=TRUE,inset=-0.085,ncol=2)
}
