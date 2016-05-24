VS.Allele.Hist <-
function(Trad_case_1,Trad_case_2,Trad_control_1,Trad_control_2,
    MC_case_1,MC_case_2,MC_control_1,MC_control_2,
    Trad_col="grey",MC_col="black",
	 main="distribution for allele frequency",title=c("case_A","case_a","control_A","control_a"),
	 xlab="count"){
	 Trad_data=rbind(Trad_case_1,Trad_case_2,Trad_control_1,Trad_control_2)
	 MC_data=rbind(MC_case_1,MC_case_2,MC_control_1,MC_control_2)
	 
    opar=par(no.readonly=TRUE)
    par(mfrow=c(2,2))
	 if(length(main)!=0){
	    par(oma=c(0, 1, 3, 0))
	 }
	 for(i in 1:4){
	    xlim_value=range(c(Trad_data[i,],MC_data[i,]))
		 Trad_ylim=max(hist(Trad_data[i,],plot=FALSE)$density)
		 MC_ylim=max(hist(MC_data[i,],plot=FALSE)$density)
		 ylim_value=c(0,max(Trad_ylim,MC_ylim))
		 
       hist(Trad_data[i,],freq=FALSE,border=Trad_col,
	       main=title[i],xlim=xlim_value,ylim=ylim_value,col=Trad_col,xlab=xlab)
	 
	    par(new=TRUE)
	    hist(MC_data[i,],freq=FALSE,border=MC_col,
	       main="",xlim=xlim_value,ylim=ylim_value,xlab="")
			 
	    legend("top",legend=c("TradPerm","MCPerm"), fill=c(Trad_col,"white"),
		    border=c("white",MC_col),cex=0.8,ncol=2,
			 xpd=TRUE,inset=-0.22,
			 box.col="white",text.font=1,merge=FALSE)
	 }
	 mtext(text=main,side=3,outer=TRUE)
	 par(opar)
}
