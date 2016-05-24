VS.Genotype.CDC <-
function(Trad_case_11,Trad_case_12,Trad_case_22,
    Trad_control_11,Trad_control_12,Trad_control_22,
	 MC_case_11,MC_case_12,MC_case_22,MC_control_11,MC_control_12,MC_control_22,
    Trad_col="black",MC_col="red",title=NULL,xlab="Genotype count",ylab="cumulative probability"){
    
	 title_value=c("case_AA","case_Aa","case_aa","control_AA","control_Aa","control_aa")
	 
	 Trad_data=rbind(Trad_case_11,Trad_case_12,Trad_case_22,
	    Trad_control_11,Trad_control_12,Trad_control_22)
	 MC_data=rbind(MC_case_11,MC_case_12,MC_case_22,
	    MC_control_11,MC_control_12,MC_control_22)
	 
	 opar=par(no.readonly=TRUE)
    par(mfrow=c(2,3))
	 if(length(title)!=0){
	    par(oma=c(0, 1, 3, 0))
	 }
	 for(i in 1:6){
	    xlim_value=range(c(Trad_data[i,],MC_data[i,]))
		 plot(ecdf(Trad_data[i,]),col=Trad_col,main=title_value[i],xlim=xlim_value,
		    xlab=xlab,ylab=ylab)	
	 
	    par(new=TRUE)
	    plot(ecdf(MC_data[i,]),col=MC_col,main='',xlim=xlim_value,
		    xlab="",ylab="")
 
	    legend("topleft",legend=c("TradPerm","MCPerm"), fill=c(Trad_col,MC_col),
		    border=c(Trad_col,MC_col))
	 }
	 mtext(text=title,side=3,outer=TRUE)
	 par(opar)
}
