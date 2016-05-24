VS.CDC <-
function(Trad_data,MC_data,Trad_col="black",MC_col="red",
    title=NULL,xlab=NULL,ylab="cumulative probability"){
		 
    xlim_value=range(c(Trad_data,MC_data))
	 plot(ecdf(Trad_data),col=Trad_col,main=title,xlim=xlim_value,
		    xlab=xlab,ylab=ylab)	
	 
	 par(new=TRUE)
	 plot(ecdf(MC_data),col=MC_col,main='',xlim=xlim_value,
		    xlab="",ylab="")
 
	 legend("top",legend=c("TradPerm","MCPerm"), fill=c(Trad_col,MC_col),
		 border=c(Trad_col,MC_col),xpd=TRUE,inset=-0.12,ncol=2,box.col="white")
}
