VS.Genotype.QQ <-
function(Trad_case_11,Trad_case_12,Trad_case_22,
    Trad_control_11,Trad_control_12,Trad_control_22,
	 MC_case_11,MC_case_12,MC_case_22,MC_control_11,MC_control_12,MC_control_22,
    scatter_col="black",line_col="black",title=NULL,
	 xlab="Quantile of genotype count (TradPerm)",
	 ylab="Quantile of genotype count (MCPerm)"){
	 
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
	    qqplot(Trad_data[i,],MC_data[i,],main=title_value[i],col=scatter_col,
	       xlab=xlab,
	       ylab=ylab)
	    abline(a=0,b=1,lwd=1,col=line_col)	 
	 }
	 mtext(text=title,side=3,outer=TRUE)
	 par(opar)
}
