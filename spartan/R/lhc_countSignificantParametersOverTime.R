lhc_countSignificantParametersOverTime<-function(FILEPATH,MEASURES,TIMEPOINTS)
{
	#Count the number of signficant parameters over time (where P Value < 0.01)
	
	for(m in 1:length(MEASURES))
	{
		SIGPARAMS_OVER_TIME<-NULL
		
		MEASURE<-MEASURES[m]
		
		PVALS<-read.csv(paste(FILEPATH,"/All_Timepoint_PVALS_",MEASURE,".csv",sep=""),header=TRUE)
		
		for(t in 1:length(TIMEPOINTS))
		{
			TIMEPOINTPROCESSING<-TIMEPOINTS[t]
			COLHEAD<-paste(MEASURE,"_PValue_",TIMEPOINTPROCESSING,sep="")
			# Read in the column for this timepoint, and see what is below 0.01
			T_PVALS<-PVALS[COLHEAD]
			
			SIGPARAMS_OVER_TIME<-rbind(SIGPARAMS_OVER_TIME,(cbind(TIMEPOINTPROCESSING,nrow(subset(T_PVALS,T_PVALS[,COLHEAD]<0.01)))))
			
		}
		
		# Plot the number of parameters that are significant for this measure
		png(filename=paste(FILEPATH,"/",MEASURE,"_Significant_Parameters.png",sep=""))
		plot(TIMEPOINTS,SIGPARAMS_OVER_TIME[,2],type="l",ylim=c(0,max(SIGPARAMS_OVER_TIME[,2])),col="red",xlab="Hours",ylab="Significant Parameters",main=MEASURE)
		abline(v=500,lty=2,col="gray")
		abline(v=1000,lty=2,col="gray")
		abline(v=1500,lty=2,col="gray")
		abline(h=5,lty=2,col="gray")
		abline(h=10,lty=2,col="gray")
		abline(h=15,lty=2,col="gray")
		abline(h=20,lty=2,col="gray")
		legend("bottomright", legend=c("Parameters with p<0.01"),lty=c(1),lwd=c(2.5),col=c("red"))
		dev.off()
	}
	
	
	
}
