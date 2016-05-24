# Plots PRCC over time for each parameter against that of the dummy parameter
lhc_graphPRCCForMultipleTimepoints<-function(FILEPATH,MEASURES,TIMEPOINTS)
{
	# For all measures
	for(m in 1:length(MEASURES))
	{
		MEASURE<-MEASURES[m]
		
		# Read in the PRCC summary
		ALL_PRCCS<- read.csv(paste(FILEPATH,"/All_Timepoint_PRCCS_",MEASURE,".csv",sep=""),header=T)
		
		# NOW TO GRAPH EACH PARAMETER(ROW) - MINUS 1 TO NOT DO THIS FOR THE DUMMY
		for(PARAM in 1:(nrow(ALL_PRCCS)-1))
		{
			# GOING TO GRAPH THIS ROW, AND THE DUMMY (WHICH SHOULD BE THE FINAL ROW)
			PARAMDATA<-ALL_PRCCS[PARAM, 2:ncol(ALL_PRCCS) ]
			DUMMYDATA<-ALL_PRCCS[nrow(ALL_PRCCS), 2:ncol(ALL_PRCCS) ]
			NAMES<-c(toString(ALL_PRCCS[PARAM,1]),"Dummy")
			
			png(filename=paste(FILEPATH,"/",MEASURE,"_",NAMES[1],".png",sep=""))
			plot(TIMEPOINTS,PARAMDATA,type="l",ylim=c(-1,1),col="red",xlab="Hours",ylab="Partial Rank Corrrelation Coefficient",main=paste("PRCC Over Time for Parameter ",ALL_PRCCS[PARAM, 1 ],"\n Measure: ",MEASURE,sep=""),yaxt="n")
			axis(side=2, at=seq(-1, 1, by=0.25))
			lines(TIMEPOINTS,DUMMYDATA,type="l",col="blue")
			
			abline(v=500,lty=2,col="gray")
			abline(v=1000,lty=2,col="gray")
			abline(v=1500,lty=2,col="gray")
			abline(h=0.25,lty=2,col="gray")
			abline(h=0.5,lty=2,col="gray")
			abline(h=0.75,lty=2,col="gray")
			abline(h=0.00,lty=2,col="gray")
			abline(h=-0.25,lty=2,col="gray")
			abline(h=-0.5,lty=2,col="gray")
			abline(h=-0.75,lty=2,col="gray")
			
			
			legend("topleft", legend=NAMES,lty=c(1,1),lwd=c(2.5,2.5),col=c("red","blue"))
			dev.off()
			
		}
	}
	
}