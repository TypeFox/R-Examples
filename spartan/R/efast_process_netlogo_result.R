efast_process_netlogo_result <-
function(FILEPATH,EFASTSAMPLE_RESULTFILENAME,PARAMETERS,NUMCURVES,NUMSAMPLES,MEASURES,RESULTFILENAME,TIMESTEP)
{
	# TO KEEP THIS IN LINE WITH TRADITIONAL SPARTAN, AND AS NETLOGO OFFERS THE CHANCE TO RUN REPEATED RUNS OF THE SAME EXPERIMENT
	# WE ARE GOING TO CREATE A MEDIAN SET OF RESULTS FOR EACH LHC SAMPLE. THEN THE OLD SPARTAN CODE CAN BE USED FROM THAT POINT
	
	# BUT FIRST, WE MAY NEED TO CHANGE THE PARAMETER AND MEASURE STRINGS. WHEN R IMPORTS THE SPREADSHEET, ANY HYPHENS OR SPACES ARE CHANGED
	# TO DOTS, AND THIS NEEDS DETECTING
	# April 2015 - do not need to do this anymore, using check.names=FALSE when reading CSVs
	#PARAMETERS<-table_header_check(PARAMETERS)
	#MEASURES<-table_header_check(MEASURES)


	for(CURVENUM in 1:NUMCURVES)
	{
		for(PARAMNUM in 1:length(PARAMETERS))
		{

			print(paste("Processing Curve: ",CURVENUM," Parameter: ",PARAMETERS[PARAMNUM],sep="")) 
			# Open the parameter file
			params<-read.csv(paste(FILEPATH,"/Curve",CURVENUM,"_",PARAMETERS[PARAMNUM],".csv",sep=""),header=TRUE,check.names=FALSE)

			CURVE_PARAM_RESULT<-NULL

			for(SAMPLE in 1:NUMSAMPLES)
			{
				if(file.exists(paste(FILEPATH,"/",CURVENUM,"/",PARAMNUM,"/",SAMPLE,"/",EFASTSAMPLE_RESULTFILENAME,SAMPLE,".csv",sep="")))
				{
					print(paste("Processing EFAST Results for Curve: ",CURVENUM," Parameter: ", PARAMNUM," Value Sample: ",SAMPLE,sep=""))

					# READ IN THE RESULT FILE
					# SKIP THE FIRST 6 LINES AS NONE OF THIS INFORMATION IS REQUIRED
					NL_RESULT<-read.csv(paste(FILEPATH,"/",CURVENUM,"/",PARAMNUM,"/",SAMPLE,"/",EFASTSAMPLE_RESULTFILENAME,SAMPLE,".csv",sep=""),sep=",",skip=6,check.names=FALSE)

					# ORDER IT BY RUN FOR EFFICIENCY LATER
					NL_RESULT_ORDERED<-NL_RESULT[order(NL_RESULT[,1]),]
	
					# REMOVE ALL THE OTHER TIMESTEPS AS NOT REQUIRED
					# THE TIMESTEP IS IN THE COLUMN HEADED X.step.
					# SET THIS TO NULL TO PLEASE CRAN SUBMISSION NOTE, THEN SUBSET
					# KA: REMOVED THIS APRIL 2015
					#X.step.<-NULL
					#TIMESTEP_RESULTS<-subset(NL_RESULT_ORDERED,X.step.==TIMESTEP)
					TIMESTEP_RESULTS<-subset(NL_RESULT_ORDERED,NL_RESULT_ORDERED["[step]"]==TIMESTEP)

					# NOW TO CREATE THE RESULTS FOR THIS SAMPLE SET
					# NETLOGO DOES GIVE THE OPTION OF RUNNING REPLICATES OF THE SAME EXPERIMENT
					# SO THERE MAY BE A FEW ROWS HERE. THE SUMMARY METHOD WILL SUMMARISE THESE
					
					# FIRST LETS SET UP THE NUMBER OF PARAMETER ROWS
					param_set<-params[SAMPLE,]

					# Make duplicates of the parameters to match the number of replicate runs	
					PARAMS<-NULL
					for(paramval in 1:ncol(param_set))
					{
						PARAMS<-cbind(PARAMS,param_set[[paramval]])
					}

					
					DUP_PARAMS<-NULL
					for(r in 1:nrow(TIMESTEP_RESULTS)-1)
					{
						DUP_PARAMS<-rbind(DUP_PARAMS,PARAMS)
					}

					# NOW WE CAN ADD THE RESULTS FOR EACH NETLOGO RUN
					for(RESPONSE in 1:length(MEASURES))
					{
						DUP_PARAMS<-cbind(DUP_PARAMS,TIMESTEP_RESULTS[MEASURES[RESPONSE]][,1])
					}

					CURVE_PARAM_RESULT<-rbind(CURVE_PARAM_RESULT,DUP_PARAMS)				
				}
				else
				{
					print(paste("ERROR: Results for Curve ",CURVENUM," Parameter ",PARAMNUM," Sample ",SAMPLE," not found",sep=""))

				}
			}

			colnames(CURVE_PARAM_RESULT)<-c(colnames(params),MEASURES)

			# Write this file out to the FILEPATH
			RESULTSFILE = paste(FILEPATH,"/Curve",CURVENUM,"_Parameter",PARAMNUM,"_Results.csv",sep="")
			write.csv(CURVE_PARAM_RESULT,RESULTSFILE,quote = FALSE,row.names=FALSE)
		}

		print(paste("Analysis of Netlogo Results for Curve ",CURVENUM," Complete",sep=""))
	}
}
	

