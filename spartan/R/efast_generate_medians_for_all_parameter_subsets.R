efast_generate_medians_for_all_parameter_subsets <-
function(FILEPATH,NUMCURVES,PARAMETERS,NUMSAMPLES,NUMRUNSPERSAMPLE,MEASURES,RESULTFILENAME,ALTERNATIVEFILENAME,OUTPUTCOLSTART,
		OUTPUTCOLEND,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL)
{
	# SPARTAN VERSION 2.0 - THIS FILE NOW MAKES A MEDIANS FILE PER PARAMETER, PER CURVE, TO FIT IN WITH THE NEW WAY THAT RESULTS CAN BE 
	# PROVIDED IN A CSV FILE. FOR AN AGENT-BASED SIMULATION, THIS WILL BE COMPRISED OF THE MEDIAN OF EACH MEASURE FOR EACH RUN

	if(is.null(TIMEPOINTS))
	{
		if(file.exists(FILEPATH))
		{
			print("Generating Simulation Median Response Sets (efast_generate_medians_for_all_parameter_subsets)")
			for(CURVE in 1:NUMCURVES)				# REPRESENTS THE CURVES
			{
				# NOW LOOK AT EACH PARAMETER OF INTEREST
				for(PARAM in 1:length(PARAMETERS))
				{
					print(paste("Generating Median Simulation Results for Curve: ",CURVE," Parameter: ",PARAM,sep=""))

					# Open the parameter file
					params<-read.csv(paste(FILEPATH,"/Curve",CURVE,"_Param",PARAM,".csv",sep=""),header=TRUE,check.names=FALSE)

					CURVE_PARAM_RESULT<-NULL

					for(j in 1:NUMSAMPLES)
					{

						SAMPLEFILEDIR<-paste(FILEPATH,"/",CURVE,"/",PARAM,"/",j,"/",sep="")
	
						medians<-getMediansSubset(SAMPLEFILEDIR,NUMRUNSPERSAMPLE,MEASURES,RESULTFILENAME,ALTERNATIVEFILENAME,OUTPUTCOLSTART,OUTPUTCOLEND)

						if(!is.null(medians))
						{
							# GET THE ROW OF PARAMETERS FROM THE FILE
							param_set<-params[j,]

							# Make duplicates of the parameters to match the number of replicate runs	
							PARAMS<-NULL
							for(paramval in 1:ncol(param_set))
							{
								PARAMS<-cbind(PARAMS,param_set[[paramval]])
							}

							DUP_PARAMS<-NULL
							for(r in 1:nrow(medians)-1)
							{
								DUP_PARAMS<-rbind(DUP_PARAMS,PARAMS)
							}

							# Now combine medians with paramters
							RESULT<-cbind(DUP_PARAMS,medians)
					
							# ADD TO THE LIST OF ALL 65 RESULTS
							CURVE_PARAM_RESULT<-rbind(CURVE_PARAM_RESULT,RESULT)
						}
						else
						{
							print(paste(CURVE," ",PARAM,sep=""))
						}
					}

					colnames(CURVE_PARAM_RESULT)<-c(colnames(params),MEASURES)

					# Write this file out to the FILEPATH
					RESULTSFILE = paste(FILEPATH,"/Curve",CURVE,"_Parameter",PARAM,"_Results.csv",sep="")
					write.csv(CURVE_PARAM_RESULT,RESULTSFILE,quote = FALSE,row.names=FALSE)
		
				}
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

			if(!is.null(ALTERNATIVEFILENAME))
			{
				ALTERNATIVEFILENAMEFULL<-paste(substr(ALTERNATIVEFILENAME,0,nchar(ALTERNATIVEFILENAME)-4),"_",TIMEPOINTPROCESSING,".",RESULTFILEFORMAT,sep="")
			}
			else
			{
				ALTERNATIVEFILENAMEFULL<-ALTERNATIVEFILENAME
			}

			efast_generate_medians_for_all_parameter_subsets(FILEPATH,NUMCURVES,PARAMETERS,NUMSAMPLES,NUMRUNSPERSAMPLE,
									MEASURES,SIMRESULTFILENAME,ALTERNATIVEFILENAMEFULL,OUTPUTCOLSTART,
									OUTPUTCOLEND,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL)
		}
	}
}

