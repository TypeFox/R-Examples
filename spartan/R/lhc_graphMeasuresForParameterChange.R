lhc_graphMeasuresForParameterChange <-
function(FILEPATH,PARAMETERS,MEASURES,MEASURE_SCALE,CORCOEFFSOUTPUTFILE,LHCSUMMARYFILENAME,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL)
{
	if(is.null(TIMEPOINTS) || length(TIMEPOINTS)==1)
	{
		if(file.exists(FILEPATH))
		{
			# LHCSUMMARYFILENAME IS LHCSummary.csv FOR 1 TIMEPOINT
			# CORCOEFFSOUTPUTFILE IS corCoefs.csv FOR 1 TIMEPOINT
			if(file.exists(paste(FILEPATH,"/",CORCOEFFSOUTPUTFILE,sep="")))
			{
				CORCOEFFS<-read.csv(paste(FILEPATH,"/",CORCOEFFSOUTPUTFILE,sep=""),header=TRUE,check.names=FALSE)
				
				# Not using this anymore, as using check.names=FALSE on CSV file reading
				# Check the Measures and Parameters for Spaces - R will have replaced these with a dot
				#MEASURES<-table_header_check(MEASURES)

				if(file.exists(paste(FILEPATH,"/",LHCSUMMARYFILENAME,sep="")))
				{
					LHCRESULTFILE<-read.csv(paste(FILEPATH,"/",LHCSUMMARYFILENAME,sep=""),header=TRUE,check.names=FALSE)

					print("Generating output graphs for LHC Parameter Analysis (lhc_graphMeasuresForParameterChange)")
	
					# CREATE A GRAPH FOR EACH PARAMETER, FOR EACH MEASURE
					for(p in 1:length(PARAMETERS))
					{
						for(m in 1:length(MEASURES))
						{
							# CREATE LABELS
							yLabel<-paste("Median Value Across Runs ",MEASURE_SCALE[m],sep="")
							xLabel<-"Parameter Value"
							# CREATE CORRELATION LABEL FOR ABOVE GRAPH
							correlationLab<-paste(MEASURES[m],"_Estimate",sep="")
							# GET THE CORRELATION FIGURE
							corrResult<-CORCOEFFS[p,correlationLab]
			
							# Where the resulting graph should go
							if(is.null(TIMEPOINTS))
							{
								GRAPHFILE <- paste(FILEPATH,"/",PARAMETERS[p],"_",MEASURES[m],".pdf",sep="")
								GRAPHTITLE <- paste("LHC Analysis for Parameter: ",PARAMETERS[p],"\nMeasure: ",MEASURES[m],
								"\nCorrelation Coefficient: ",toString(signif(corrResult,3)),sep="")
							}
							else
							{
								GRAPHFILE <- paste(FILEPATH,"/",PARAMETERS[p],"_",MEASURES[m],"_",TIMEPOINTS,".pdf",sep="")
								GRAPHTITLE <- paste("LHC Analysis for Parameter: ",PARAMETERS[p],"\nMeasure: ",MEASURES[m],
								"\nCorrelation Coefficient: ",toString(signif(corrResult,3))," Timepoint: ",TIMEPOINTS," ",TIMEPOINTSCALE,sep="")
							}
			
							pdf(GRAPHFILE)
			
							plot(LHCRESULTFILE[,PARAMETERS[p]],LHCRESULTFILE[,MEASURES[m]],type="p",pch=4,xlab=xLabel,ylab=yLabel,
								main=GRAPHTITLE)
						
							# NO LEGEND NECESSARY AS STATED IN THE HEADER
							par(xpd=FALSE)
			
							dev.off()
						}
					}
					print("LHC Graphs Complete")		
				}
				else
				{
					print("Cannot find LHC Summary File. Are you sure you have run the method to generate it?")
				}
			}
			else
			{
				print("Cannot find Partial Rank Correlation Coefficients File. Are you sure you have run the method to generate it?")
			}
		}
		else
		{
			print("The directory specified in FILEPATH does not exist. No Output Graphs Generated")
		}
	}
	else
	{
		# PROCESS EACH TIMEPOINT, BY AMENDING THE FILENAMES AND RECALLING THIS FUNCTION
		for(n in 1:length(TIMEPOINTS))
		{
			TIMEPOINTPROCESSING<-TIMEPOINTS[n]
			print(paste("PROCESSING TIMEPOINT: ",TIMEPOINTPROCESSING,sep=""))

			CORCOEFFSOUTPUTFILE_FORMAT<-substr(CORCOEFFSOUTPUTFILE,(nchar(CORCOEFFSOUTPUTFILE)+1)-3,nchar(CORCOEFFSOUTPUTFILE))
			CORCOEFFSOUTPUTFILE_FULL<-paste(substr(CORCOEFFSOUTPUTFILE,0,nchar(CORCOEFFSOUTPUTFILE)-4),"_",TIMEPOINTPROCESSING,".",CORCOEFFSOUTPUTFILE_FORMAT,sep="")	

			LHCSUMMARYFILENAME_FORMAT<-substr(LHCSUMMARYFILENAME,(nchar(LHCSUMMARYFILENAME)+1)-3,nchar(LHCSUMMARYFILENAME))
			LHCSUMMARYFILENAME_FULL<-paste(substr(LHCSUMMARYFILENAME,0,nchar(LHCSUMMARYFILENAME)-4),"_",TIMEPOINTPROCESSING,".",LHCSUMMARYFILENAME_FORMAT,sep="")

			lhc_graphMeasuresForParameterChange(FILEPATH,PARAMETERS,MEASURES,MEASURE_SCALE,CORCOEFFSOUTPUTFILE_FULL,LHCSUMMARYFILENAME_FULL,
								TIMEPOINTS=TIMEPOINTPROCESSING,TIMEPOINTSCALE=TIMEPOINTSCALE)

		}

	}
}

