VS.Allele.CDC <-
function(Trad_case_1,Trad_case_2,Trad_control_1,Trad_control_2,
    MC_case_1,MC_case_2,MC_control_1,MC_control_2,
    Trad_col="black",MC_col="red",main="cumulative distribution curve",
	 title=c("case_A","case_a","control_A","control_a"),
	 xlab="count",ylab="cumulative probability"){
	 
	 Trad_data=rbind(Trad_case_1,Trad_case_2,Trad_control_1,Trad_control_2)
	 MC_data=rbind(MC_case_1,MC_case_2,MC_control_1,MC_control_2)
	 
	 opar=par(no.readonly=TRUE)
    par(mfrow=c(2,2))
	 if(length(title)!=0){
	    par(oma=c(0, 1, 3, 0))
	 }
	 for(i in 1:4){
	    xlim_value=range(c(Trad_data[i,],MC_data[i,]))
		 plot(ecdf(Trad_data[i,]),col=Trad_col,main=title[i],xlim=xlim_value,
		    xlab=xlab,ylab=ylab)	
	 
	    par(new=TRUE)
	    plot(ecdf(MC_data[i,]),col=MC_col,main='',xlim=xlim_value,
		    xlab="",ylab="")
 
	    legend("topleft",legend=c("TradPerm","MCPerm"), fill=c(Trad_col,MC_col),
		    border=c(Trad_col,MC_col),cex=0.7)
	 }
	 mtext(text=main,side=3,outer=TRUE)
	 par(opar)
}
