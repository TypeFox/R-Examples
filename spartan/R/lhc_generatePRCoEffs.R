lhc_generatePRCoEffs <-
function(FILEPATH,PARAMETERS,MEASURES,LHCSUMMARYFILENAME,CORCOEFFSOUTPUTFILE,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL)
{
	if(is.null(TIMEPOINTS))
	{
		if(file.exists(FILEPATH))
		{
			# LHCSUMMARYFILENAME IS LHCSummary.csv FOR 1 TIMEPOINT
			# CORCOEFFSOUTPUTFILE IS corCoefs.csv FOR 1 TIMEPOINT
			if(file.exists(paste(FILEPATH,"/",LHCSUMMARYFILENAME,sep="")))
			{
				LHCRESULTFILE<-read.csv(paste(FILEPATH,"/",LHCSUMMARYFILENAME,sep=""),header=TRUE, check.names=FALSE)
				PARAMCOEFFSTRUCT<-NULL
				COEFFRESULTS<-NULL
				print("Generating Partial Rank Correlation Coefficients (lhc_generatePRCoEffs)")
				
				# Not using anymore as we're using check.names=FALSE on CSV file reading
				# Check the Measures and Parameters for Spaces - R will have replaced these with a dot
				#MEASURES<-table_header_check(MEASURES)
		
				# NEED TO GENERATE A COEFFICIENT FOR EACH PARAMETER BEING EXAMINED
				for(k in 1:length(PARAMETERS))
				{
					PARAMNAME<-PARAMETERS[k]
		
					# GET COEFFICIENT SET
					COEFFDATA<-lhc_constructCoEffDataSet(LHCRESULTFILE,PARAMNAME,PARAMETERS)
					# GET PARAMETER RESULT
					COEFFPARAMCOL<-LHCRESULTFILE[,PARAMETERS[k]]
		
					PARAMRESULTS<-NULL
					# GET MEASURE RESULTS AND CALCULATE COEFFICIENTS FOR EACH PARAMETER
					for(l in 1:length(MEASURES))
					{
						COEFFMEASURERESULT<-LHCRESULTFILE[,MEASURES[l]]
						PARAMCOEFF<-pcor.test(COEFFPARAMCOL,COEFFMEASURERESULT,COEFFDATA,calcMethod=c("s"))
						PARAMRESULTS<-cbind(PARAMRESULTS,PARAMCOEFF$estimate,PARAMCOEFF$p.value)
					}
		
					COEFFRESULTS<-rbind(COEFFRESULTS,PARAMRESULTS)		
				}
			
				# NAME THE COLUMNS FOR EASE OF REFERENCE LATER
				COEFFRESULTSHEAD<-NULL	
				for(l in 1:length(MEASURES))
				{
					COEFFRESULTSHEAD<-cbind(COEFFRESULTSHEAD,
								(paste(MEASURES[l],"_Estimate",sep="")),
								(paste(MEASURES[l],"_PValue",sep="")))
				}
	
				# OUTPUT THE RESULTS FOR ALL PARAMETERS
				colnames(COEFFRESULTS)<-c(COEFFRESULTSHEAD)
				rownames(COEFFRESULTS)<-PARAMETERS
		
				COEFFSRESULTSFILE<-paste(FILEPATH,"/",CORCOEFFSOUTPUTFILE,sep="")
				write.csv(COEFFRESULTS,COEFFSRESULTSFILE,quote = FALSE)
		
				print(paste("File of Partial Rank Correlation Coefficients Generated. Output to ",COEFFSRESULTSFILE,sep=""))
			}
			else
			{
				print("LHC Summary file cannot be found. Are you sure you have run the method to generate it?")
			}
		}
		else
		{
			print("The directory specified in FILEPATH does not exist. No Partial Rank Correlation Coefficients Generated")
		}
	}
	else
	{
		# PROCESS EACH TIMEPOINT, BY AMENDING THE FILENAMES AND RECALLING THIS FUNCTION
		for(n in 1:length(TIMEPOINTS))
		{
			TIMEPOINTPROCESSING<-TIMEPOINTS[n]
			print(paste("Processing Timepoint: ",TIMEPOINTPROCESSING,sep=""))

			LHCSUMMARYFILENAME_FORMAT<-substr(LHCSUMMARYFILENAME,(nchar(LHCSUMMARYFILENAME)+1)-3,nchar(LHCSUMMARYFILENAME))
			LHCSUMMARYFILENAME_FULL<-paste(substr(LHCSUMMARYFILENAME,0,nchar(LHCSUMMARYFILENAME)-4),"_",TIMEPOINTPROCESSING,".",LHCSUMMARYFILENAME_FORMAT,sep="")

			CORCOEFFSOUTPUTFILE_FORMAT<-substr(CORCOEFFSOUTPUTFILE,(nchar(CORCOEFFSOUTPUTFILE)+1)-3,nchar(CORCOEFFSOUTPUTFILE))
			CORCOEFFSOUTPUTFILE_FULL<-paste(substr(CORCOEFFSOUTPUTFILE,0,nchar(CORCOEFFSOUTPUTFILE)-4),"_",TIMEPOINTPROCESSING,".",CORCOEFFSOUTPUTFILE_FORMAT,sep="")

			lhc_generatePRCoEffs(FILEPATH,PARAMETERS,MEASURES,LHCSUMMARYFILENAME_FULL,CORCOEFFSOUTPUTFILE_FULL,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL)
		}
	}

}

