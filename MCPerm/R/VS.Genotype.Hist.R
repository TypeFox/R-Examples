VS.Genotype.Hist <-
function(Trad_case_11,Trad_case_12,Trad_case_22,
    Trad_control_11,Trad_control_12,Trad_control_22,
	 MC_case_11,MC_case_12,MC_case_22,MC_control_11,MC_control_12,MC_control_22,
    Trad_col="grey",MC_col="black",title=NULL,xlab="Genotype count"){
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
		 Trad_ylim=max(hist(Trad_data[i,],plot=FALSE)$density)
		 MC_ylim=max(hist(MC_data[i,],plot=FALSE)$density)
		 ylim_value=c(0,max(Trad_ylim,MC_ylim))
		 
       hist(Trad_data[i,],freq=FALSE,border=Trad_col,
	       main=title_value[i],xlim=xlim_value,ylim=ylim_value,col=Trad_col,xlab=xlab)
	 
	    par(new=TRUE)
	    hist(MC_data[i,],freq=FALSE,border=MC_col,
	       main="",xlim=xlim_value,ylim=ylim_value,xlab="")
			 
	    legend("topright",legend=c("TradPerm","MCPerm"), fill=c(Trad_col,"white"),
		    border=c("white",MC_col),cex=0.7)
	 }
	 mtext(text=title,side=3,outer=TRUE)
	 par(opar)
}
