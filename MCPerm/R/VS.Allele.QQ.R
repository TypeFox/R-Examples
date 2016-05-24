VS.Allele.QQ <-
function(Trad_case_1,Trad_case_2,Trad_control_1,Trad_control_2,
    MC_case_1,MC_case_2,MC_control_1,MC_control_2,
    scatter_col="black",line_col="black",
	 main="QQ plot for allele model",
	 title=c("case_A","case_a","control_A","control_a"),
	 xlab="Quantile of count (TradPerm)",
	 ylab="Quantile of count (MCPerm)"){
	 
	 Trad_data=rbind(Trad_case_1,Trad_case_2,Trad_control_1,Trad_control_2)
	 MC_data=rbind(MC_case_1,MC_case_2,MC_control_1,MC_control_2)
	 
	 opar=par(no.readonly=TRUE)
    par(mfrow=c(2,2))
	 if(length(title)!=0){
	    par(oma=c(0, 1, 3, 0))
	 }
	 for(i in 1:4){
	    qqplot(Trad_data[i,],MC_data[i,],main=title[i],col=scatter_col,
	       xlab=xlab,
	       ylab=ylab)
	    abline(a=0,b=1,lwd=1,col=line_col)	 
	 }
	 mtext(text=main,side=3,outer=TRUE)
	 par(opar)
}
