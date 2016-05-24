lhc_calculatePRCCForMultipleTimepoints<-function(FILEPATH, CORCOEFFSOUTPUTFILE, TIMEPOINTS, MEASURES)
{
	# Calculates the PRCC for each parameter at each timepoint
	# Unlike former Spartan, this stores PRCC and P-Value in 2 different files to make the plot function easier 
	
	for(m in 1:length(MEASURES))
	{
		MEASURE<-MEASURES[m]
		PRCC_LABEL<-paste(MEASURE,"_Estimate",sep="")
		PVAL_LABEL<-paste(MEASURE,"_PValue",sep="")
		
		PRCCS_OVER_TIME<-NULL
		PVALS_OVER_TIME<-NULL
		PRCC_HEADERS<-c("Parameter Name")
		PVALS_HEADERS<-c("Parameter Name")
		
		for(t in 1:length(TIMEPOINTS))
		{
			TIMEPOINTPROCESSING<-TIMEPOINTS[t]
			CORCOEFFSOUTPUTFILE_FORMAT<-substr(CORCOEFFSOUTPUTFILE,(nchar(CORCOEFFSOUTPUTFILE)+1)-3,nchar(CORCOEFFSOUTPUTFILE))
			CORCOEFFSOUTPUTFILE_FULL<-paste(substr(CORCOEFFSOUTPUTFILE,0,nchar(CORCOEFFSOUTPUTFILE)-4),"_",TIMEPOINTPROCESSING,".",CORCOEFFSOUTPUTFILE_FORMAT,sep="")
			
			if(file.exists(paste(FILEPATH,"/",CORCOEFFSOUTPUTFILE_FULL,sep="")))
			{
				COEFFS_TIMEPOINT<-read.csv(paste(FILEPATH,"/",CORCOEFFSOUTPUTFILE_FULL,sep=""),header=T)
								
				if(t==1)
				{
					# Copy over the parameter name in this instance as well as the results
					PRCCS_OVER_TIME<-COEFFS_TIMEPOINT["X"]
					PRCCS_OVER_TIME<-cbind(PRCCS_OVER_TIME,COEFFS_TIMEPOINT[PRCC_LABEL])
					PVALS_OVER_TIME<-COEFFS_TIMEPOINT["X"]
					PVALS_OVER_TIME<-cbind(PVALS_OVER_TIME,COEFFS_TIMEPOINT[PVAL_LABEL])
				}else
				{
					PRCCS_OVER_TIME<-cbind(PRCCS_OVER_TIME,COEFFS_TIMEPOINT[PRCC_LABEL])
					PVALS_OVER_TIME<-cbind(PVALS_OVER_TIME,COEFFS_TIMEPOINT[PVAL_LABEL])
				}
				
				PRCC_HEADERS<-cbind(PRCC_HEADERS,paste(PRCC_LABEL,"_",TIMEPOINTPROCESSING,sep=""))
				PVALS_HEADERS<-cbind(PVALS_HEADERS,paste(PVAL_LABEL,"_",TIMEPOINTPROCESSING,sep=""))
			}
			else
			{
				print(paste("Correlation Coefficients file for Timepoint ",TIMEPOINTPROCESSING," does not exist",sep=""))
			}
			
		}
		
		# Write the summary to file, if not empty
		if(!is.null(PRCCS_OVER_TIME))
		{
			# ADD HEADERS TO THE PRCC RESULTS
			colnames(PRCCS_OVER_TIME)<-PRCC_HEADERS
			colnames(PVALS_OVER_TIME)<-PVALS_HEADERS
			
			RESULTSFILE = paste(FILEPATH,"/All_Timepoint_PRCCS_",MEASURE,".csv",sep="")
			write.csv(PRCCS_OVER_TIME,RESULTSFILE,quote = FALSE,row.names=FALSE)
			
			RESULTSFILE = paste(FILEPATH,"/All_Timepoint_PVALS_",MEASURE,".csv",sep="")
			write.csv(PVALS_OVER_TIME,RESULTSFILE,quote = FALSE,row.names=FALSE)
		}
		else
		{
			print("No Correlation Coefficients to write to file")
		}
	}
		
}