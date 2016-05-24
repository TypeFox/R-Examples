plotATestsFromTimepointFiles<-function(FILEPATH,PARAMETERS,ATESTRESULTFILENAME,ATESTSIGLEVEL,MEASURES,PMIN,PMAX,PINC,TIMEPOINTS)
{
	for(PARAM in 1:length(PARAMETERS))
	{
		FULLPARAMRESULTS<-data.frame()
		CNAME<-NULL
		
		for(i in 1:length(TIMEPOINTS))
		{
			hour<-TIMEPOINTS[i]
			
			OATResults<-read.csv(paste(FILEPATH,PARAMETERS[PARAM],"/",ATESTRESULTFILENAME,"_",hour,".csv",sep=""),header=T)
			
			# ADD THE PARAMETER VALUES
			if(ncol(FULLPARAMRESULTS)==0)
			{
				FULLPARAMRESULTS<-OATResults[1]
				CNAME<-c("ParameterValue")
			}
			
			# Add the result for all measures
			for(m in 1:length(MEASURES))
			{
				FULLPARAMRESULTS<-cbind(FULLPARAMRESULTS,OATResults[paste("ATest",MEASURES[m],sep="")])
				CNAME<-cbind(CNAME,paste(MEASURES[m],"_",hour,sep=""))
			}
		}
		
		colnames(FULLPARAMRESULTS)<-CNAME
		
		# Now produce the A-Tests timepoints plot for all measures
		
		for(j in 1:length(MEASURES))
		{
			# PLOT EACH TIMEPOINT
			# DO THE FIRST
			GRAPHFILE <- paste(FILEPATH,"/",PARAMETERS[PARAM],"_",MEASURES[j],".pdf",sep="")
			pdf(GRAPHFILE,width=7,height=7.2)
			
			MEASURELABEL<-paste(MEASURES[j],"_",TIMEPOINTS[1],sep="")
			
			GRAPHTITLE<-paste("One-A-Time Parameter Analysis Over Simulation Time\nParameter: ",
					PARAMETERS[PARAM],", Measure: ",MEASURES[j],sep="")
			
			plot(FULLPARAMRESULTS$ParameterValue,FULLPARAMRESULTS[,MEASURELABEL],type="o",main=GRAPHTITLE,lty=1,ylim=c(0,1),pch=1,xlab = "Parameter Value",ylab = "A Test Score",xaxt="n")
			
			# NOW ADD THE REST
			
			for(l in 2:length(TIMEPOINTS))
			{
				MEASURELABEL<-paste(MEASURES[j],"_",TIMEPOINTS[l],sep="")
				lines(FULLPARAMRESULTS$ParameterValue,FULLPARAMRESULTS[,MEASURELABEL],type="o",lty=5,pch=l)
			}
			
			# Add the x axis
			axis(1,at=seq(PMIN[PARAM],PMAX[PARAM],by=PINC[PARAM]))
			
			# legend	
			legend("topleft", inset=.0,title="Timepoints",TIMEPOINTS, pch=1:length(TIMEPOINTS),cex=0.75)
			
			par(xpd=FALSE)
			
			# Lines for A-Test scores
			abline(a=0.5,b=0,lty=4)
			text((PMAX[PARAM]+PMIN[PARAM])/2, 0.52, "no difference", col = "blue") 
			abline(a=(0.5+ATESTSIGLEVEL),b=0,lty=4)
			text((PMAX[PARAM]+PMIN[PARAM])/2,(0.5+ATESTSIGLEVEL+0.02), "large difference", col = "blue") 
			abline(a=(0.5-ATESTSIGLEVEL),b=0,lty=4)
			text((PMAX[PARAM]+PMIN[PARAM])/2,(0.5-ATESTSIGLEVEL-0.02), "large difference", col = "blue") 
			dev.off()
		}
	}
}
