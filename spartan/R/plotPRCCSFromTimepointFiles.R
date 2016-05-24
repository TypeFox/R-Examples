plotPRCCSFromTimepointFiles<-function(FILEPATH,PARAMETERS,MEASURES,CORCOEFFSFILENAME,TIMEPOINTS,TIMEPOINTSCALE,DISPLAYPVALS=FALSE)
{
	print("Plotting Graphs for Partial Rank Correlation Coefficients Over Time")
	
	if(requireNamespace("plotrix",quietly=TRUE))
	{
		# One plot for each parameter
		for(PARAM in 1:length(PARAMETERS))
		{
			# PRCCS for this parameter
			FULLPARAMRESULTS<-NULL
			# P-Values for this parameter
			PARAMPVALS<-data.frame()
			
			# Now to gather the data for each hour from the relevant result files
			for(i in 1:length(TIMEPOINTS))
			{
				hour<-TIMEPOINTS[i]
				
				# Add the timepoint to the correlation coefficient results file
				CORCOEFFSOUTPUTFILE_FORMAT<-substr(CORCOEFFSFILENAME,(nchar(CORCOEFFSFILENAME)+1)-3,nchar(CORCOEFFSFILENAME))
				CORCOEFFSOUTPUTFILE_FULL<-paste(substr(CORCOEFFSFILENAME,0,nchar(CORCOEFFSFILENAME)-4),"_",hour,".",CORCOEFFSOUTPUTFILE_FORMAT,sep="")
				
				
				# Read in the coefficients
				LHCResults<-read.csv(paste(FILEPATH,"/",CORCOEFFSOUTPUTFILE_FULL,sep=""),header=T)
				# Get the PRCCS
				results<-c(hour,LHCResults[PARAM,2],LHCResults[PARAM,4])		
				# Get the P-Values
				pvals.d<-data.frame(LHCResults[PARAM,3],LHCResults[PARAM,5])
				
				# Append the PRCCS for this timepoint to those of all timepoints
				FULLPARAMRESULTS<-rbind(FULLPARAMRESULTS,results)
				# Append the P-Values for this timepoint to those of all timepoints
				PARAMPVALS<-rbind(PARAMPVALS,pvals.d)
			}
			
			# Set the row and column names of the P-Values data frame
			rownames(PARAMPVALS)<-TIMEPOINTS
			colnames(PARAMPVALS)<-MEASURES
			
			# Now to make the plot
			GRAPHFILE <- paste(FILEPATH,"/",PARAMETERS[PARAM],"_OverTime.pdf",sep="")
			pdf(GRAPHFILE,width=7,height=7)
			
			# Title, with parameter name
			GRAPHTITLE<-paste("Partial Rank Correlation Coefficients Over Simulation Time\nParameter: ",
					PARAMETERS[PARAM],sep="")
			
			# Plot the first measure	
			plot(FULLPARAMRESULTS[,1],FULLPARAMRESULTS[,2],type="o",main=GRAPHTITLE,lty=1,xlab="",ylim=c(-1,1),ylab = "Partial Rank Correlation Coefficient",xaxt="n",yaxt="n",bty="n")
			
			# Now add the rest
			lines(FULLPARAMRESULTS[,1],FULLPARAMRESULTS[,3],type="o",lty=5,pch=2)
			
			axis(2,at=seq(-1,1,by=0.25))
			axis(1,pos=0,at=seq(as.numeric(min(TIMEPOINTS)),as.numeric(max(TIMEPOINTS)),by=(as.numeric(max(TIMEPOINTS))/length(TIMEPOINTS))),tck=0.015,labels=FALSE)
			# Add the axis at 0
			abline(h=0)
			# Add the labels to the axis
			for(h in 1:length(TIMEPOINTS))
			{
				text(as.numeric(TIMEPOINTS[h]),0.08,TIMEPOINTS[h])
			}
			
			# Add the X axis label
			text(((as.numeric(max(TIMEPOINTS))/2)+as.numeric(min(TIMEPOINTS))/2),0.18,TIMEPOINTSCALE)
			
			# P-Values Table, if the user wants this displayed
			if(DISPLAYPVALS==TRUE)
			{
				xAxisLoc=(((as.numeric(max(TIMEPOINTS))-as.numeric(min(TIMEPOINTS)))/100)*71)+as.numeric(min(TIMEPOINTS))
				plotrix::addtable2plot(xAxisLoc,0.7,signif(PARAMPVALS,digits=3),cex=0.7,display.rownames=TRUE,title="p-Values",display.colnames=TRUE,bty="o",hlines=TRUE)
			}
			
			# Graph legend
			legend("topleft", inset=.025,title="Measures",MEASURES, pch=1:length(MEASURES))
			
			# Output graph
			dev.off()
			
		}
		print(paste("Complete. Check for output in the directory ",FILEPATH,sep=""))
	}
	else
	{
		print("The plotPRCCSFromTimepointFiles function requires the plotrix package to be installed")
	}
}
