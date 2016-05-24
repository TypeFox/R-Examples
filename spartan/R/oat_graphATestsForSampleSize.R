oat_graphATestsForSampleSize <-
function(FILEPATH,PARAMETERS,MEASURES,ATESTSIGLEVEL,ATESTRESULTFILENAME,BASELINE,PMIN=NULL,PMAX=NULL,PINC=NULL,PARAMVALS=NULL,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL)
{
	if(is.null(TIMEPOINTS) || length(TIMEPOINTS)==1)
	{
		# NOTE THAT OUTPUT_FOLDER AND BASELINE PARAMETERS ADDED IN SPARTAN 2.0
		# IN VERSION 2.0, WE PROCESS ONE FILE, NOT AN ATEST FILE FOR EACH PARAMETER
		# THE GRAPHS ALSO GO IN THE SAME DIRECTORY AS THIS FOLDER, NOT IN THE TOP LEVEL

		if(file.exists(FILEPATH))
		{
			print("Creating graphs of A-Test results (oat_graphATestsForSampleSize)") 

			# FIRSTLY READ IN THE ATESTS FILE (NOW ALL IN ONE FILE)
			RESULT<-read.csv(paste(FILEPATH,"/",ATESTRESULTFILENAME,sep=""),sep=",",header=TRUE, check.names=FALSE)

			# No longer in use as have changed the method of reading in CSV files (check.names=FALSE)
			# Check the Measures and Parameters for Spaces - R will have replaced these with a dot
			#MEASURES<-table_header_check(MEASURES)

			for(PARAM in 1:length(PARAMETERS))
			{
				print(paste("Creating graph for Parameter ",PARAMETERS[PARAM],sep=""))

				#if(file.exists(paste(FILEPATH,"/",PARAMETERS[PARAM],"/",ATESTRESULTFILENAME,".csv",sep="")))
				#{
				# NOW WE NEED TO RECOVER THE A-TESTS FOR THIS PARAMETER FROM THE SET OF RESULTS
				EXP_PARAMS<-as.character(BASELINE)

				# NOW GET THE LIST OF PARAMETER VALUES BEING EXPLORED FOR THIS PARAMETER
				# NOTE THE CONVERSION BACK TO NUMBERS - GETS RID OF TRAILING ZEROS MADE BY SEQ
				PARAM_VAL_LIST<-as.numeric(prepare_parameter_value_list(PMIN,PMAX,PINC,PARAMVALS,PARAM))

				PARAM_ATESTS<-NULL

				# NOW WE ITERATE THROUGH THE VALUES IN THIS LIST
				for(PARAMVAL in 1:length(PARAM_VAL_LIST))
				{
					# NOW TO RECOVER THE EXPERIMENTS RUN UNDER THIS PARAMETER SET
					ATESTS<-RESULT

					# SET THE VALUE OF THIS PARAMETER TO BE THAT WE ARE PROCESSING
					EXP_PARAMS[PARAM] <- as.character(PARAM_VAL_LIST[PARAMVAL])

					for(PARAMOFINT in 1:length(PARAMETERS))
					{
						ATESTS<-subset(ATESTS,ATESTS[[PARAMETERS[PARAMOFINT]]]==as.numeric(EXP_PARAMS[PARAMOFINT]))
					}

					# KEEP ALL THE ATESTS RESULTS TOGETHER FOR THIS PARAMETER
					PARAM_ATESTS<-rbind(PARAM_ATESTS,ATESTS)

				}
	
				# Where the resulting graph should go 
				if(is.null(TIMEPOINTS))
				{
					GRAPHFILE <- paste(FILEPATH,"/",PARAMETERS[PARAM],".pdf",sep="")
					GRAPHTITLE <- paste("A-Test Scores when adjusting parameter \n",PARAMETERS[PARAM],sep="")
				}
				else
				{
					GRAPHFILE <- paste(FILEPATH,"/",PARAMETERS[PARAM],"_",TIMEPOINTS,".pdf",sep="")
					GRAPHTITLE <- paste("A-Test Scores when adjusting parameter \n",PARAMETERS[PARAM]," at Timepoint: ",TIMEPOINTS," ",TIMEPOINTSCALE,sep="")
				}
				
				#pdf(GRAPHFILE,width=12,height=7)
				#par(xpd=NA,oma=c(0,0,0,14))
				pdf(GRAPHFILE, width=12, height=7)
				par(xpd=NA,mar=c(4,4,4,17))
		
				# NOW PLOT THE MEASURES
				# START WITH THE FIRST
				MEASURELABEL<-paste("ATest",MEASURES[1],sep="")
				plot(PARAM_ATESTS[[PARAMETERS[PARAM]]],PARAM_ATESTS[,MEASURELABEL],type="o",main=GRAPHTITLE,lty=1,ylim=c(0,1),pch=1,xlab = "Parameter Value",ylab = "A Test Score",xaxt="n")
		
				if(length(MEASURES)>1)
				{
					# NOW ADD THE REST OF THE MEASURES
					for(l in 2:length(MEASURES))
					{
						MEASURELABEL<-paste("ATest",MEASURES[l],sep="")
						#lines(ATESTS$ParameterVal,ATESTS[,MEASURELABEL],type="o",lty=5,pch=l)
						lines(PARAM_ATESTS[[PARAMETERS[PARAM]]],PARAM_ATESTS[,MEASURELABEL],type="o",lty=5,pch=l)
					}
				}
			
				axis(1,PARAM_VAL_LIST)
				#legend(par("usr")[2],par("usr")[4],title="Measures",MEASURES,pch=1:length(MEASURES),lty=1,xjust=0,yjust=2.0)
				legend(par("usr")[2],par("usr")[4],title="Measures",MEASURES, pch=1:length(MEASURES),cex=0.7,ncol=1)
				par(xpd=FALSE)
		
				abline(a=0.5,b=0,lty=4)
				text((max(PARAM_VAL_LIST)+min(PARAM_VAL_LIST))/2, 0.52, "no difference", col = "blue") 
				abline(a=(0.5+ATESTSIGLEVEL),b=0,lty=4)
				text((max(PARAM_VAL_LIST)+min(PARAM_VAL_LIST))/2,(0.5+ATESTSIGLEVEL+0.02), "large difference", col = "blue") 
				abline(a=(0.5-ATESTSIGLEVEL),b=0,lty=4)
				text((max(PARAM_VAL_LIST)+min(PARAM_VAL_LIST))/2,(0.5-ATESTSIGLEVEL-0.02), "large difference", col = "blue") 
	
				dev.off()
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

			ATESTRESULTFILENAME_FORMAT<-substr(ATESTRESULTFILENAME,(nchar(ATESTRESULTFILENAME)+1)-3,nchar(ATESTRESULTFILENAME))
			ATESTRESULTFILENAME_FULL<-paste(substr(ATESTRESULTFILENAME,0,nchar(ATESTRESULTFILENAME)-4),"_",TIMEPOINTPROCESSING,".",ATESTRESULTFILENAME_FORMAT,sep="")
			
			oat_graphATestsForSampleSize(FILEPATH,PARAMETERS,MEASURES,ATESTSIGLEVEL,ATESTRESULTFILENAME_FULL,BASELINE,
							PMIN,PMAX,PINC,PARAMVALS,TIMEPOINTS=TIMEPOINTPROCESSING,TIMEPOINTSCALE=TIMEPOINTSCALE)
		}
	}
}

