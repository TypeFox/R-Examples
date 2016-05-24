oat_plotResultDistribution <-
function(FILEPATH,PARAMETERS,MEASURES,MEASURE_SCALE,CSV_FILE_NAME,BASELINE,PMIN=NULL,PMAX=NULL,PINC=NULL,PARAMVALS=NULL,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL)
{
	if(is.null(TIMEPOINTS) || length(TIMEPOINTS)==1)
	{
		if(file.exists(FILEPATH))
		{
			print("Plotting result distribution for each parameter (oat_plotResultDistribution)")

			# NEW TO SPARTAN VERSION 2
			# READS SIMULATION RESPONSES FROM A CSV FILE, IN THE FORMAT: PARAMETER VALUES (COLUMNS), SIMULATION OUTPUT MEASURES
			# IN A CHANGE TO SPARTAN 1, THE FIRST FUNCTION THAT PROCESSES SIMULATION RESPONSES CREATES THIS FILE, NOT MEDIANS FOR EACH PARAMETER AS IT USED TO 
			# THIS WAY WE ARE NOT DEALING WITH TWO METHODS OF SIMULATION RESULT SPECIFICATION
			# READ IN THE OAT RESULT FILE
			RESULT<-read.csv(paste(FILEPATH,"/",CSV_FILE_NAME,sep=""),sep=",",header=TRUE, check.names=FALSE)
			
			# No longer in use as have changed the method of reading CSV files (check.names=FALSE)
			# Check the Measures and Parameters for Spaces - R will have replaced these with a dot
			#MEASURES<-table_header_check(MEASURES)

			for(PARAM in 1:length(PARAMETERS))
			{
				print(paste("Creating Output Responses Box Plot Graph for Parameter ",PARAMETERS[PARAM],sep=""))

				# THE RESULTS OF THE OAT ANALYSIS IS IN ONE PLACE. THUS WE NEED TO REFER TO THE CORRECT BASELINE RESULT FOR PARAMETERS THAT ARE
				# NOT BEING CHANGED SO WE USE THE VARIABLE EXP_PARAMS WHEN WE START A NEW VARIABLE - WE SET THE PARAMS TO THE BASELINE AND THEN ONLY
				# ALTER THE ONE BEING CHANGED
				EXP_PARAMS<-as.character(BASELINE)

				# NOW GET THE LIST OF PARAMETER VALUES BEING EXPLORED FOR THIS PARAMETER
				# NOTE THE CONVERSION BACK TO NUMBERS - GETS RID OF TRAILING ZEROS MADE BY SEQ
				PARAM_VAL_LIST<-as.numeric(prepare_parameter_value_list(PMIN,PMAX,PINC,PARAMVALS,PARAM))

				ALLRESULTS<-NULL

				for(PARAMVAL in 1:length(PARAM_VAL_LIST))
				{
					EXP_PARAMS[PARAM] <- as.character(PARAM_VAL_LIST[PARAMVAL])
					PARAM_RESULT<-subset_results_by_param_value_set(PARAMETERS,RESULT,EXP_PARAMS)

					VALUE_RESULT<-cbind(PARAM_RESULT[PARAMETERS[PARAM]])
					
					# NOW ADD ALL MEASURES
					for(MEASURE in 1:length(MEASURES))
					{
						VALUE_RESULT<-cbind(VALUE_RESULT,PARAM_RESULT[MEASURES[MEASURE]])
					}

					ALLRESULTS<-rbind(ALLRESULTS,VALUE_RESULT)
				}

				for(MEASURE in 1:length(MEASURES))
				{
					# NOW DO THE BOXPLOTS FOR EACH MEASURE
					# BOXPLOT THE MEASURE
					if(is.null(TIMEPOINTS))
					{
						GRAPHFILE = paste(FILEPATH,"/",PARAMETERS[PARAM],"_",
							MEASURES[MEASURE],"_BP.pdf",sep="")
						GRAPHTITLE<-paste("Distribution of ",MEASURES[MEASURE], " Responses \n when altering parameter ",
							PARAMETERS[PARAM],sep="")
					}else
					{
						GRAPHFILE = paste(FILEPATH,"/",PARAMETERS[PARAM],"_",
							MEASURES[MEASURE],"_BP_",TIMEPOINTS,".pdf",sep="")
						GRAPHTITLE<-paste("Distribution of ",MEASURES[MEASURE], " Responses \n when altering parameter ",
							PARAMETERS[PARAM]," at Timepoint ",TIMEPOINTS," ",TIMEPOINTSCALE,sep="")
					}
		
					pdf(GRAPHFILE)			
		
					# GENERATE YLABEL BASED ON PARAMETER MEASURE
					YLABEL<-paste("Median ",MEASURES[MEASURE]," (",MEASURE_SCALE[PARAM],")",sep="")
					MEASURESLAB<-MEASURES[MEASURE]
				
					boxplot(ALLRESULTS[,MEASURE+1]~ALLRESULTS[,1],ylab=YLABEL,xlab="Parameter Value",main=GRAPHTITLE)
				
					dev.off()

					print(paste("Box Plot Generated and output as ",GRAPHFILE,sep=""))

				}
			}
		}
		else
		{
			print("The directory specified in FILEPATH does not exist. No graph created")
		}
	}
	else
	{
		# PROCESS EACH TIMEPOINT, BY AMENDING THE FILENAMES AND RECALLING THIS FUNCTION
		for(n in 1:length(TIMEPOINTS))
		{
			TIMEPOINTPROCESSING<-TIMEPOINTS[n]
			print(paste("PROCESSING TIMEPOINT: ",TIMEPOINTPROCESSING,sep=""))

			CSV_FILE_NAME_FORMAT<-substr(CSV_FILE_NAME,(nchar(CSV_FILE_NAME)+1)-3,nchar(CSV_FILE_NAME))
			CSV_FILE_NAME_FULL<-paste(substr(CSV_FILE_NAME,0,nchar(CSV_FILE_NAME)-4),"_",TIMEPOINTPROCESSING,".",CSV_FILE_NAME_FORMAT,sep="")

			oat_plotResultDistribution(FILEPATH,PARAMETERS,MEASURES,MEASURE_SCALE,CSV_FILE_NAME_FULL,BASELINE,PMIN,
							PMAX,PINC,PARAMVALS,TIMEPOINTS=TIMEPOINTPROCESSING,TIMEPOINTSCALE=TIMEPOINTSCALE)

		}
	}

}

