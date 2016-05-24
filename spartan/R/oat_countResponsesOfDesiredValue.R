oat_countResponsesOfDesiredValue <-
function(FILEPATH,PARAMETERS,RESULTFILENAME,OUTPUTCOLSTART,OUTPUTCOLEND,PARAMETER,NUMRUNSPERSAMPLE,MEASURE,DESIREDRESULT,OUTPUTFILENAME,BASELINE,PMIN=NULL,PMAX=NULL,PINC=NULL,PARAMVALS=NULL,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL)
{
	if(is.null(TIMEPOINTS))
	{
		if(file.exists(FILEPATH))
		{
			print(paste("Summing Responses for Parameter ",PARAMETER," where ", MEASURE, " = ",DESIREDRESULT,sep=""))

			if(file.exists(paste(FILEPATH,"/",PARAMETER,sep="")))
			{

				EXP_PARAMS<-as.character(BASELINE)

				# NOW WE CAN WORK WITH INCREMENTS BETWEEN MAX AND MIN, AND SPECIFIED VALUES, WE NEED TO GET THE VALUES
				# OF THE PARAMETERS WE ARE ANALYSING
				# NOTE THE CONVERSION BACK TO NUMBERS - GETS RID OF TRAILING ZEROS MADE BY SEQ
				PARAM_VAL_LIST<-as.numeric(prepare_parameter_value_list(PMIN,PMAX,PINC,PARAMVALS,match(PARAMETER,PARAMETERS)))

				ALLRESULTS<-NULL

				# NOW WE ITERATE THROUGH THE VALUES IN THIS LIST
				for(PARAMVAL in 1:length(PARAM_VAL_LIST))
				{
					# SET THE VALUE OF THE PARAMETERS BEING EXAMINED TO INCLUDE THE CURRENT VALUE OF THE PARAMETER
					EXP_PARAMS[match(PARAMETER,PARAMETERS)] <- as.character(PARAM_VAL_LIST[PARAMVAL])

					TRUECOUNT=0
					FALSECOUNT=0
			
					if(file.exists(paste(FILEPATH,"/",PARAMETER,"/",toString(PARAMVAL),sep="")))
					{
						for(i in 1:NUMRUNSPERSAMPLE)
						{
							FILEADDRESS = paste(FILEPATH,"/",PARAMETER,"/",toString(PARAMVAL),"/",i,"/",RESULTFILENAME,sep="")
	
							if(RESULTFILEFORMAT=="csv")
							{
								# DEALING WITH A CSV FILE
								if(file.exists(paste(FILEADDRESS,".",RESULTFILEFORMAT,sep="")))
								{
									if(OUTPUTCOLSTART>1)
									{
										import<-read.csv(paste(FILEADDRESS,".csv",sep=""),colClasses=c(rep('NULL',OUTPUTCOLSTART-1),rep(NA,OUTPUTCOLEND-OUTPUTCOLSTART+1)), header=TRUE, check.names=FALSE)
									}else
									{
										import<-read.csv(paste(FILEADDRESS,".csv",sep=""),colClasses=c(rep(NA,OUTPUTCOLEND)), header=TRUE, check.names=FALSE)
									}

									MODELRESULT<-data.frame(import,check.names=FALSE)

								}else if(RESULTFILEFORMAT=="xml")
								{
									if(requireNamespace("XML",quietly=TRUE))
									{
										MODELRESULT<-XML::xmlToDataFrame(paste(FILEADDRESS,".xml",sep=""))
									}
									else
									{
										print("The oat_countResponsesOfDesiredValue function requires the XML package to be installed")
									}
								}


								if(MODELRESULT[[MEASURE]]==DESIREDRESULT)
								{
									TRUECOUNT=TRUECOUNT+1
								}else
								{
									FALSECOUNT=FALSECOUNT+1
								}								
							}
						}
					
						ROWRESULT<-cbind(PARAMVAL,TRUECOUNT,FALSECOUNT)
						ALLRESULTS<-rbind(ALLRESULTS,ROWRESULT)
					}
					else
					{
						print(paste("No results can be found for parameter: ",PARAMETER," Value: ",PARAMVAL,sep=""))
					}
				}

				# Output the true false results
				RESULTSFILE = paste(FILEPATH,"/",PARAMETER,"/",MEASURE,"_",DESIREDRESULT,sep="")
				write.csv(ALLRESULTS,paste(RESULTSFILE,".csv",sep=""),quote = FALSE,row.names=FALSE)
			}
			else
			{
				print(paste("No results can be found for the parameter specified: ",PARAMETER,sep=""))
			}
		
		}
		else
		{
			print("The directory specified in FILEPATH does not exist. No analysis completed")
		}
	}
	else
	{
		# PROCESS EACH TIMEPOINT, BY AMENDING THE FILENAMES AND RECALLING THIS FUNCTION
		for(n in 1:length(TIMEPOINTS))
		{

			TIMEPOINTPROCESSING<-TIMEPOINTS[n]
			print(paste("PROCESSING TIMEPOINT: ",TIMEPOINTPROCESSING,sep=""))

			RESULTFILEFORMAT<-substr(RESULTFILENAME,(nchar(RESULTFILENAME)+1)-3,nchar(RESULTFILENAME))
			SIMRESULTFILENAME<-paste(substr(RESULTFILENAME,0,nchar(RESULTFILENAME)-4),"_",TIMEPOINTPROCESSING,".",RESULTFILEFORMAT,sep="")
		
			OUTPUTFILENAMEFORMAT<-substr(OUTPUTFILENAME,(nchar(OUTPUTFILENAME)+1)-3,nchar(OUTPUTFILENAME))
			OUTPUTFILENAME_FULL<-paste(substr(OUTPUTFILENAME,0,nchar(OUTPUTFILENAME)-4),"_",TIMEPOINTPROCESSING,".",OUTPUTFILENAMEFORMAT,sep="")

			# NOW CALL THIS FUNCTION AGAIN TO DO THE TIMEPOINTS - WE SET THE TIMEPOINTS AND TIMEPOINTSCALE TO NULL NOW SO WE DONT END UP BACK IN THIS ELSE
			oat_countResponsesOfDesiredValue(FILEPATH,SIMRESULTFILENAME,OUTPUTCOLSTART,OUTPUTCOLEND,PARAMETER,NUMRUNSPERSAMPLE,MEASURE,DESIREDRESULT,OUTPUTFILENAME_FULL,
							BASELINE,PMIN,PMAX,PINC,PARAMVALS,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL)
		}
	}

}
