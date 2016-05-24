lhc_plotCoEfficients<-function(FILEPATH, CORCOEFFSOUTPUTFILE, MEASURES, PRINTOPT, TIMEPOINTS=NULL, TIMEPOINTSCALE=NULL)
{
	if(is.null(TIMEPOINTS) || length(TIMEPOINTS)==1)
	{
		if(file.exists(FILEPATH))
		{
			# READ IN THE COEFFICIENTS FILE
			COEFFSRESULTSFILENAME<-paste(FILEPATH,"/",CORCOEFFSOUTPUTFILE,sep="")
			COEFFS<-read.csv(COEFFSRESULTSFILENAME,header=TRUE,check.names=FALSE)
			
			# Check the Measures and Parameters for Spaces - R will have replaced these with a dot
			#MEASURES<-table_header_check(MEASURES)
			
			# COLUMN 1 HAS PARAMETER NAMES, THEN FOLLOWS FOR EACH MEASURE - THE PRCC AND THE P VALUE
			# WE'RE GOING TO GRAPH ALL THE PRCC'S ON ONE GRAPH
			
			if(PRINTOPT=="INDIVIDUAL")
			{
				# INDIVIDUAL PLOTS FOR EACH MEASURE
				print("Producing Partial Rank Correlation Coefficient Plots for Each Measure")
				
				for(i in 1:length(MEASURES))
				{	
					if(is.null(TIMEPOINTS))
					{
						GRAPHFILE <- paste(FILEPATH,"/PRCC_Measure_",MEASURES[i],".pdf",sep="")
						GRAPHTITLE <- paste("PRCC Values for Measure: ",MEASURES[i],sep="")
					}
					else
					{
						GRAPHFILE <- paste(FILEPATH,"/PRCC_Measure_",MEASURES[i],"_",TIMEPOINTS,".pdf",sep="")
						GRAPHTITLE <- paste("PRCC Values for Measure: ",MEASURES[i],"\nTimepoint: ",TIMEPOINTS,sep="")
					}
					
					
					pdf(GRAPHFILE, width=9, height=5)
					par(xpd=NA,mar=c(2,4,2,17))
					
					# Generate the heading of the CSV file - the measure plus _Estimate
					M<-paste(MEASURES[i],"_Estimate",sep="")
					# We can now use this to get the column out the dataset
					
					barplot(COEFFS[,M],ylim=c(-1,1),col="black",
							main=GRAPHTITLE,
							ylab="Partial Rank Correlation Coefficient",
							names.arg=seq(1,nrow(COEFFS),by=1))
					
					thelabels <- paste(1:nrow(COEFFS), " ", COEFFS[,1], sep="")				
					par(xpd=TRUE)
					legend((nrow(COEFFS)+6),1.0,legend = thelabels, pch="",cex=0.7,ncol=1)
					par(xpd=FALSE)
					dev.off()
		
					
				}
				
			} else if(PRINTOPT=="ALL")
			{
				print("Producing Partial Rank Correlation Coefficient Summary Plot of All Measures")
				
				# ALL PRCCS FOR ALL MEASURES, ON ONE PLOT
				# Make the data frame for the plot
				# FIRST OF ALL WE NEED TO REMOVE THE P VALUES SO WE CAN AUTOMATE THIS
				
				if(is.null(TIMEPOINTS))
				{
					GRAPHFILE <- paste(FILEPATH,"/PRCC_AllMeasures.pdf",sep="")
					GRAPHTITLE <- "PRCC Values for All Measures"
				}
				else
				{
					GRAPHFILE <- paste(FILEPATH,"/PRCC_AllMeasures_",TIMEPOINTS,".pdf",sep="")
					GRAPHTITLE <- paste("PRCC Values for All Measures\nTimepoint: ",TIMEPOINTS,sep="")
				}
				
				pdf(GRAPHFILE,width=9, height=5)
				
				par(xpd=NA,mar=c(2,4,2,17))
				
				PRCCS<-NULL
				for(p in seq(2,ncol(COEFFS),by=2))
				{
					PRCCS<-cbind(PRCCS,COEFFS[,p])
				}
				
				# NOW MAKE THE DATA FRAME
				d<-data.frame(row.names=levels(COEFFS[,1]),PRCCS,check.names=FALSE)
				colnames(d)<-MEASURES
				d <- do.call(rbind, d)
				barplot(d, beside = TRUE, ylim=c(-1,1.4), legend.text = rownames(d),
						args.legend = list(x = "topright", bty="n"),names.arg=seq(1,nrow(COEFFS),by=1),
						main=GRAPHTITLE,ylab="Partial Rank Correlation Coefficient")
				
				thelabels <- paste(1:nrow(COEFFS), " ", COEFFS[,1], sep="")				
				par(xpd=TRUE)
				legend((nrow(COEFFS)+6),1.0,legend = thelabels, pch="",cex=0.7,ncol=1)
				par(xpd=FALSE)
				dev.off()
				
				#leg<-NULL
				#for(v in seq(1,nrow(COEFFS),by=1))
				#{
				#	leg<-cbind(leg,toString(v))
				#}
				
				#par(xpd=TRUE)
				#legend(0.9,-1.25, COEFFS[,1],pch=leg)
				
				#par(xpd=FALSE)
				
				#dev.off()
				
			}
		}
		else
		{
			print("The directory specified in FILEPATH does not exist. No analysis performed")
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
			
			lhc_plotCoEfficients(FILEPATH, CORCOEFFSOUTPUTFILE_FULL, MEASURES, PRINTOPT, TIMEPOINTPROCESSING,NULL)
			
		}
	}
}
