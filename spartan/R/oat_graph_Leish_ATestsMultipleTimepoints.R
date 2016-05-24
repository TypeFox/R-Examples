oat_graph_Leish_ATestsMultipleTimepoints <-
		function(FILEPATH,PARAMETERS,MEASURES,PMIN,PMAX,PINC,PARAMVALS,BASELINE,ATESTRESULTFILENAME,ATESTSIGLEVEL,TIMEPOINTS)
{
	for(PARAM in 1:length(PARAMETERS))
	{
		print(paste("Creating graph for Parameter ",PARAMETERS[PARAM],sep=""))
		
		for(m in 1:length(MEASURES))
		{
			MEASURE<-MEASURES[m]
			MEASURE_ATESTS<-NULL
			
			for(t in 1:length(TIMEPOINTS))
			{
				# Read in the timepoint file
				
				ATESTRESULTFILENAME_FORMAT<-substr(ATESTRESULTFILENAME,(nchar(ATESTRESULTFILENAME)+1)-3,nchar(ATESTRESULTFILENAME))
				ATESTRESULTFILENAME_FULL<-paste(substr(ATESTRESULTFILENAME,0,nchar(ATESTRESULTFILENAME)-4),"_",TIMEPOINTS[t],".",ATESTRESULTFILENAME_FORMAT,sep="")
				
				RESULT<-read.csv(paste(FILEPATH,ATESTRESULTFILENAME_FULL,sep=""),header=T)
				
				#print(t)
				#if(file.exists(paste(FILEPATH,"/",PARAMETERS[PARAM],"/",ATESTRESULTFILENAME,".csv",sep="")))
				#{
				# NOW WE NEED TO RECOVER THE A-TESTS FOR THIS PARAMETER FROM THE SET OF RESULTS
				EXP_PARAMS<-as.character(BASELINE)
				
				# NOW GET THE LIST OF PARAMETER VALUES BEING EXPLORED FOR THIS PARAMETER
				# NOTE THE CONVERSION BACK TO NUMBERS - GETS RID OF TRAILING ZEROS MADE BY SEQ
				PARAM_VAL_LIST<-as.numeric(prepare_parameter_value_list(PMIN,PMAX,PINC,PARAMVALS,PARAM))
				
				if(t==1)
				{
					MEASURE_ATESTS<-c(PARAM_VAL_LIST)
				}
				
				PARAM_ATESTS<-NULL
				
				# NOW WE ITERATE THROUGH THE VALUES IN THIS LIST
				for(PARAMVAL in 1:length(PARAM_VAL_LIST))
				{
					#print(PARAMVAL)
					# NOW TO RECOVER THE EXPERIMENTS RUN UNDER THIS PARAMETER SET
					ATESTS<-RESULT
					
					# SET THE VALUE OF THIS PARAMETER TO BE THAT WE ARE PROCESSING
					EXP_PARAMS[PARAM] <- as.character(PARAM_VAL_LIST[PARAMVAL])
					
					for(PARAMOFINT in 1:length(PARAMETERS))
					{
						ATESTS<-subset(ATESTS,ATESTS[[PARAMETERS[PARAMOFINT]]]==as.numeric(EXP_PARAMS[PARAMOFINT]))
						#print(nrow(ATESTS))
					}
					
					# KEEP ALL THE ATESTS RESULTS TOGETHER FOR THIS PARAMETER
					PARAM_ATESTS<-rbind(PARAM_ATESTS,ATESTS[1,])
					
				}	
				
				# ADDED: ONLY KEEP WHAT WE REALLY NEED FOR THIS ANALYSIS
				LABEL<-paste("ATest",MEASURE,sep="")
				LABEL2<-paste("ATest",MEASURE,"Norm",sep="")
				ATESTS<-cbind(PARAM_ATESTS[LABEL],PARAM_ATESTS[LABEL2])
				colnames(ATESTS)<-c(paste(MEASURE,"_",TIMEPOINTS[t],sep=""),paste(MEASURE,"_Norm_",TIMEPOINTS[t],sep=""))
				MEASURE_ATESTS<-cbind(MEASURE_ATESTS,ATESTS)
			}
			
			# Now plot this parameter/measure pair
			# Start with the first timepoint
			#GRAPHFILE <- paste(FILEPATH,"/",PARAMETERS[PARAM],"_",MEASURE,".pdf",sep="")
			png(filename=paste(FILEPATH,"/",PARAMETERS[PARAM],"_",MEASURE,".png",sep=""))
			LABEL<-paste(MEASURE,"_",TIMEPOINTS[1],sep="")	
			#pdf(GRAPHFILE, width=12, height=7)
			par(xpd=NA,mar=c(4,4,4,4))
			
			GRAPHTITLE <- paste("A-Test Scores when adjusting parameter ",PARAMETERS[PARAM],"\n Measure: ",MEASURE,sep="")
			plot(MEASURE_ATESTS[,1],MEASURE_ATESTS[[LABEL]],type="o",main=GRAPHTITLE,lty=1,ylim=c(0,1),pch=1,xlab = "Parameter Value",ylab = "A Test Score",xaxt="n")
			
			# Now add all the other timepoints
			for(t in 2:length(TIMEPOINTS))
			{
				LABEL<-paste(MEASURE,"_",TIMEPOINTS[t],sep="")	
				lines(MEASURE_ATESTS[,1],MEASURE_ATESTS[[LABEL]],type="o",lty=1,pch=t)
			}
			
			axis(1,PARAM_VAL_LIST)
			legend(par("usr")[2],par("usr")[4],title="Timepoints",TIMEPOINTS, pch=1:length(TIMEPOINTS),cex=0.7,ncol=1)
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
}