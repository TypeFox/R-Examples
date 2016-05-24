lhc_generate_lhc_sample_netlogo <-
function(FILEPATH,PARAMETERS,PARAMVALS,NUMSAMPLES,ALGORITHM,EXPERIMENT_REPETITIONS,RUNMETRICS_EVERYSTEP,NETLOGO_SETUP_FUNCTION,NETLOGO_RUN_FUNCTION,MEASURES)
{
	if(requireNamespace("lhs",quietly=TRUE))
	{
		# CHECK THE OUTPUT FILE PATH EXISTS
		if(file.exists(FILEPATH))
		{
			# AS WITH TRADITIONAL SPARTAN, WE ARE ALSO GOING TO MAKE A SPREADSHEET CONTAINING A SUMMARY OF THE PARAMETERS FOR EACH RUN
			# AS WELL AS NETLOGO EXPERIMENT FILES. INITIALISE THE STORES HERE
			LHC_PARAMETERS_ALL_RUNS<-NULL
			LHC_CSV_HEADERS<-NULL

			# Firstly we need to calculate the number of parameters that are being varied
			# We do this by detecting the number that are in square brackets and have min & max separated by commas
			NUM_PARAMSOFINT<-0
			for(PARAM in 1:length(PARAMETERS))
			{
			

				PARAMVALSPLIT<-strsplit(PARAMVALS[PARAM],",")[[1]]
	
				if(length(PARAMVALSPLIT)>1)
				{
					NUM_PARAMSOFINT <- NUM_PARAMSOFINT+1 

					# ADD THE PARAMETER TO THE CSV FILE HEADER
					LHC_CSV_HEADERS<-cbind(LHC_CSV_HEADERS,PARAMETERS[PARAM])
				}

			
			}

			# NOW WE CAN GENERATE A LHC SAMPLE FOR THIS NUMBER OF PARAMETERS
			# FIRSTLY MAKE SURE THE ALGORITHM CHOICE IS SPECIFIED IN LOWER CASE	
			ALGORITHM<-tolower(ALGORITHM)
			LHC_DESIGN<-NULL

			# PERFORM THE SAMPLING - JUDGING ON USERS CHOICE OF ALGORITHM
			if(ALGORITHM=="optimum")
			{
				# MAY TAKE A WHILE FOR A LARGE NUMBER OF SAMPLES (THIS BEING TWO DAYS WHERE NUMSAMPLES=500)
				LHC_DESIGN<-lhs::optimumLHS(NUMSAMPLES,NUM_PARAMSOFINT,maxSweeps=2,eps=.1)
			}else{
				LHC_DESIGN<-lhs::randomLHS(NUMSAMPLES,NUM_PARAMSOFINT)
			}

			# NOW TO CONSTRUCT THE EXPERIMENT FILE FOR EACH SET OF PARAMETERS. THUS THERE WILL BE THE SAME NUMBER OF EXPERIMENT XML
			# FILES AS SPECIFIED IN THE NUMBER OF SAMPLES

			# First check the package exists
			if(requireNamespace("XML",quietly=TRUE))
			{
				for(SAMPLE in 1:NUMSAMPLES)
				{
					# CREATE A FOLDER FOR THIS SAMPLE
					dir.create(file.path(FILEPATH,SAMPLE), showWarnings = FALSE)

					# CREATE A NEW ROW FOR THE CSV SUMMARY FILE
					LHC_SAMPLE_RUN<-NULL

					# AS SOME PARAMETERS VARY AND SOME DON'T, WE NEED TO KEEP A REFERENCE TO THE COLUMN OF THE GENERATED LHC FILE
					LHC_COL_REF<-1

					# INITIALISE THE XML FILE
					xml<-XML::xmlOutputDOM(tag="experiments")

					# NEXT TAG IN IS EXPERIMENT
					xml$addTag("experiment", attrs=c(name=paste("LHC_Sample",SAMPLE,sep=""),repetitions=EXPERIMENT_REPETITIONS, runMetricsEveryStep=RUNMETRICS_EVERYSTEP),close=FALSE)

					# NOW THE PROCEDURES TO CALL SETUP, GO, AND OUTPUT MEASURES TO ANALYSE
					xml$addTag("setup",NETLOGO_SETUP_FUNCTION)
					xml$addTag("go",NETLOGO_RUN_FUNCTION)
			
					for(MEASURE in 1:length(MEASURES))
					{
						xml$addTag("metric",MEASURES[MEASURE])	
					}


					# NOW TO SET THE VALUE OF EACH PARAMETER IN THIS RUN. SOME ARE STATIC, SOME WILL BE TAKEN FROM THE LHC
					for(PARAM in 1:length(PARAMETERS))
					{
						# SPLIT THE VALUES OF THIS PARAMETER BY THE COMMA. IF NO COMMA, THIS IS NOT A PARAMETER OF INTEREST
						# AND IS OF STATIC VALUE
						PARAMVALSPLIT<-strsplit(PARAMVALS[PARAM],",")[[1]]
	
						if(length(PARAMVALSPLIT)==1)
						{	
							# JUST GET THE VALUE - WE ADD THIS TO THE XML LATER
							VALUE=PARAMVALS[PARAM]
						}
						else
						{
							# THIS IS MORE DIFFICULT, AS WE NEED TO SET THE VALUE BASED ON THE LATIN HYPERCUBE SAMPLE
							# THE LHC GENERATES A NUMBER BETWEEN 0 AND 1, SO WE NEED TO SCALE THIS BASED ON THE MIN AND
							# MAX SPECIFIED
							# FIRST GET THE MAX AND MIN
							MIN<-as.numeric(substring(PARAMVALSPLIT[[1]],2))
							MAX<-as.numeric(substring(PARAMVALSPLIT[[2]], 1, nchar(PARAMVALSPLIT[[2]])-1))
	
							# NOW CALCULATE THE VALUE TO USE FOR THIS PARAMETER
							VALUE<-(LHC_DESIGN[SAMPLE,LHC_COL_REF] * (MAX - MIN) ) + MIN

							# INCREMENT THE LHC_COL_REF SO A DIFFERENT VARYING PARAMETER IS USED NEXT TIME
							LHC_COL_REF<-LHC_COL_REF+1

							# ADD TO THE SUMMARY
							LHC_SAMPLE_RUN<-cbind(LHC_SAMPLE_RUN,VALUE)
						}

						# NOW CREATE THE XML FOR THIS PARAMETER
						xml$addTag("enumeratedValueSet", attrs=c(variable=PARAMETERS[PARAM]),close=FALSE)
						# NOW ADD THE VALUE
						xml$addTag("value", attrs=c(value=VALUE))
	
						# CLOSE THE ENUMERATED VALUE SET TAG
						xml$closeTag()

					}

					# CLOSE THE EXPERIMENT TAG
					xml$closeTag()	

					# CLOSE THE EXPERIMENTS TAG
					xml$closeTag()	
			
					# SAVE THE XML FILE IN THE FOLDER FOR THIS EXPERIMENT
					XML::saveXML(xml,file=paste(FILEPATH,"/",SAMPLE,"/lhc_analysis_set",SAMPLE,".xml",sep=""),indent=TRUE, 
							prefix = '<?xml version="1.0" encoding="us-ascii"?>\n',
							doctype = '<!DOCTYPE experiments SYSTEM "behaviorspace.dtd">')

					# ADD THIS SAMPLE TO THE CSV SUMMARY FILE
					LHC_PARAMETERS_ALL_RUNS<-rbind(LHC_PARAMETERS_ALL_RUNS,LHC_SAMPLE_RUN)
				}
			}
			else
			{
				print("The lhc_generate_lhc_sample_netlogo function requires the XML package to be installed")
			}


			# WITH ALL SAMPLES DONE, CREATE THE CSV FILE
			colnames(LHC_PARAMETERS_ALL_RUNS)<-c(LHC_CSV_HEADERS)
			# OUTPUT THE LHC DESIGN AS A CSV FILE
			write.csv(LHC_PARAMETERS_ALL_RUNS,paste(FILEPATH,"/LHC_Parameters_for_Runs.csv",sep=""),row.names=FALSE,quote = FALSE)

			print(paste("Parameter Set Generated and Output to ",FILEPATH,"/LHC_Parameters_for_Runs.csv",sep=""))
			print(paste(NUMSAMPLES," Netlogo Experiment Files Output to ",FILEPATH,sep=""))
		}
		else
		{
			print("The directory specified in FILEPATH does not exist. No parameter samples generated")
		}	
	}
	else
	{
		print("The lhc_generate_lhc_sample_netlogo function requires the lhs package to be installed")
	}
}

