aa_graphATestsForSampleSize <-
function(FILEPATH,ATESTS,MEASURES,LARGEDIFFINDICATOR,GRAPHOUTPUTNAME,TIMEPOINT,TIMEPOINTSCALE)
{
	# Where the A-Test Results are - at present spartan only outputs test results as CSV
	#ATESTS <- read.csv(paste(FILEPATH,"/",SAMPLESIZE,"/",ATESTRESULTSFILENAME,".csv",sep=""),header=TRUE)
	ATESTS<-data.frame(ATESTS,check.names = FALSE)

	# Where the resulting graph should go
	GRAPHFILE = paste(FILEPATH,"/",GRAPHOUTPUTNAME,sep="")
	pdf(GRAPHFILE, width=12, height=7)
	par(xpd=NA,mar=c(4,4,4,17))

	# WILL PLOT EACH MEASURE IN TURN.  BUT PLOT THE INITIAL MEASURE FIRST
	MEASURELABEL<-paste("ATest",MEASURES[1],sep="")

	plot(ATESTS["Sample"][,1],ATESTS[MEASURELABEL][,1],type="o",lty=1,ylim=c(0,1),pch=1,xlab = "Run Subset / Parameter Value (Dummy)",ylab = "A Test Score",xaxt="n",xlim=c(2,20))

	# NOW DO THE REST OF THE VALUES, IF THERE IS MORE THAN ONE MEASURE
	if(length(MEASURES)>1)
	{
		for(l in 2:length(MEASURES))
		{
			MEASURELABEL<-paste("ATest",MEASURES[l],sep="")
			lines(ATESTS["Sample"][,1],ATESTS[MEASURELABEL][,1],type="o",lty=5,pch=l)	
		}
	}

	# NOW COMPLETE GRAPH
	# DETERMINE IF THIS IS BEING DONE FOR ONE TIMEPOINT OR MANY
	if(is.null(TIMEPOINT))
	{
		title(main=paste("A-Test Scores for ",nrow(ATESTS)," Dummy Parameters where \n Sample Size = ",ATESTS[1,1],sep=""))
	}
	else
	{
		# CREATES A LABEL WHICH SHOWS THE TIMEPOINT ANALYSED
		# SCALE HOLDS THE TIMEPOINT MEASURE - I.E ARE WE LOOKING AT HOURS/MINS/SECONDS/STEPS?
		title(main=paste("A-Test Scores for ",nrow(ATESTS)," Dummy Parameters where \n Sample Size = ",ATESTS[1,1]," at Timepoint: ",TIMEPOINT," ",TIMEPOINTSCALE,sep=""))
	}

	axis(1,at=seq(2,nrow(ATESTS)+1,by=2))
	legend(par("usr")[2],par("usr")[4],title="MEASURES",MEASURES, pch=1:length(MEASURES),cex=0.7,ncol=1)
	par(xpd=FALSE)

	# ADD THE SIGNIFICANCE LINES
	# FIRSTLY DOWN THE MIDDLE
	abline(a=0.5,b=0,lty=4)
	text(13, 0.52, "no difference", col = "blue") 
	# NOW ADD DIFFERENCES AS DICTATED BY USER INPUT
	abline(a=(0.5+LARGEDIFFINDICATOR),b=0,lty=4)
	text(13, (0.5+LARGEDIFFINDICATOR+0.02), "large difference", col = "blue") 
	abline(a=(0.5-LARGEDIFFINDICATOR),b=0,lty=4)
	text(13, (0.5-LARGEDIFFINDICATOR-0.02), "large difference", col = "blue") 

	dev.off()
}

