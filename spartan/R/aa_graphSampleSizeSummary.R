aa_graphSampleSizeSummary <-
function(FILEPATH,MEASURES,MAXSAMPLESIZE,SMALL,MEDIUM,LARGE,SUMMARYFILENAME,GRAPHOUTPUTFILE,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL,GRAPHLABEL=NULL)
{
	if(is.null(TIMEPOINTS))
	{
		# NOW DRAW THE GRAPH
		print("Creating Summary Graph")
		
		# No longer required as using check.names=FALSE in CSV file reading
		# Check the Measures and Parameters for Spaces - R will have replaced these with a dot
		#MEASURES<-table_header_check(MEASURES)
	
		if(file.exists(paste(FILEPATH,"/",SUMMARYFILENAME,sep="")))
		{
			aTestResults <- read.csv(paste(FILEPATH,"/",SUMMARYFILENAME,sep=""),header=TRUE,check.names=FALSE)

			# Where the resulting graph should go (ATESTMAXES.PDF USED IF ONE TIMEPOINT)
			graphFile = paste(FILEPATH,"/",GRAPHOUTPUTFILE,sep="")			
			pdf(graphFile, width=12, height=7)
			par(xpd=NA,mar=c(4,4,2,17))

			# NOW PLOT FOR EACH MEASURE
			# THE PLOT BEGINS WITH THE FIRST MEASURE
			measureLabel<-paste(MEASURES[1],"MaxA",sep="")
			plot(aTestResults$SampleSize,aTestResults[measureLabel][,1],type="o",lty=1,ylim=c(0.5,1.0),pch=1,xlab = "Sample Size",ylab = "A Test Score",xaxt="n",yaxt="n")
	
			# NOW DO ALL OTHER MEASURES, IF THERE ARE MORE THAN ONE
			if(length(MEASURES)>1)
			{
				for(l in 2:length(MEASURES))
				{
					measureLabel<-paste(MEASURES[l],"MaxA",sep="")
					lines(aTestResults$SampleSize,aTestResults[measureLabel][,1],type="o",lty=5,pch=l)
				}
			}
	
			# NOW COMPLETE GRAPH - TITLE DEPENDING ON WHETHER THIS IS ONE TIMEPOINT OR MANY
			if(is.null(GRAPHLABEL))
			{
				title("Maximum A-Test Scores for each Sample Size")
			}
			else
			{
				title(paste("Maximum A-Test Scores for each Sample Size \n Timepoint: ",GRAPHLABEL,sep=""))
			}
		
			axis(1,at=seq(0,MAXSAMPLESIZE,by=100))
			axis(2, at=seq(0.5,1.0, by=0.05))
			#legend(par("usr")[2],par("usr")[4],title="MEASURES",MEASURES,pch=1:length(MEASURES),lty=1,xjust=0,yjust=2.0)
			
			par(xpd=TRUE)
			legend(par("usr")[2],par("usr")[4],title="MEASURES",MEASURES, pch=1:length(MEASURES),cex=0.7,ncol=1)
			
			par(xpd=FALSE)
	
			# ADD THE LINES TO SHOW WHERE THE A-TEST EFFECTS ARE
			abline(h=SMALL,lty=4)
			text(MAXSAMPLESIZE/2,SMALL-0.01, "SMALL effect", col = "blue") 
			abline(h=LARGE,lty=4)
			text(MAXSAMPLESIZE/2, LARGE+0.01, "LARGE effect", col = "blue") 
			abline(h=MEDIUM,lty=4)
			text(MAXSAMPLESIZE/2,MEDIUM+0.02, "MEDIUM effect", col = "blue") 
			dev.off()

			print(paste("Summary Graph output to ",FILEPATH,"/",GRAPHOUTPUTFILE,sep=""))
		}
		else
		{
			print("Cannot find the summary file specified in parameter SUMMARYFILENAME. Check you have run the analysis to generate this file")
		}
	}
	else
	{
		for(n in 1:length(TIMEPOINTS))
		{

			TIMEPOINTPROCESSING<-TIMEPOINTS[n]
			print(paste("PROCESSING TIMEPOINT: ",TIMEPOINTPROCESSING,sep=""))

			SUMMARYFILENAME_FORMAT<-substr(SUMMARYFILENAME,(nchar(SUMMARYFILENAME)+1)-3,nchar(SUMMARYFILENAME))
			SUMMARYFILENAME_FULL<-paste(substr(SUMMARYFILENAME,0,nchar(SUMMARYFILENAME)-4),"_",TIMEPOINTPROCESSING,".",SUMMARYFILENAME_FORMAT,sep="")

			GRAPHFILENAME_FORMAT<-substr(GRAPHOUTPUTFILE,(nchar(GRAPHOUTPUTFILE)+1)-3,nchar(GRAPHOUTPUTFILE))
			GRAPHFILENAME_FULL<-paste(substr(GRAPHOUTPUTFILE,0,nchar(GRAPHOUTPUTFILE)-4),"_",TIMEPOINTPROCESSING,".",GRAPHFILENAME_FORMAT,sep="")

			GRAPHLABEL<-paste(TIMEPOINTPROCESSING," (",TIMEPOINTSCALE,")",sep="")

			aa_graphSampleSizeSummary(FILEPATH,MEASURES,MAXSAMPLESIZE,SMALL,MEDIUM,LARGE,SUMMARYFILENAME_FULL,GRAPHFILENAME_FULL,TIMEPOINTS=NULL,TIMEPOINTSCALE=NULL,GRAPHLABEL)


		}

	}		
}

