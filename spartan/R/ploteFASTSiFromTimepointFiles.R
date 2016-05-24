ploteFASTSiFromTimepointFiles<-function(FILEPATH,PARAMETERS,MEASURES,EFASTRESULTFILENAME,TIMEPOINTS,TIMEPOINTSCALE)
{
	for(m in 1:length(MEASURES))
	{
		MEASURE<-MEASURES[m]
		# Add Si onto the measure to get this from the result set
		MEASURELABEL<-paste(MEASURE,"_Si",sep="")
		
		SiMEASURERESULTS<-data.frame()
		
		for(i in 1:length(TIMEPOINTS))
		{
			hour<-TIMEPOINTS[i]
			
			# Add the timepoint onto the end of the filename
			EFASTRESULTFILENAME_FORMAT<-substr(EFASTRESULTFILENAME,(nchar(EFASTRESULTFILENAME)+1)-3,nchar(EFASTRESULTFILENAME))
			EFASTRESULTFILENAME_FULL<-paste(substr(EFASTRESULTFILENAME,0,nchar(EFASTRESULTFILENAME)-4),"_",hour,".",EFASTRESULTFILENAME_FORMAT,sep="")
			
			
			# READ IN THE TIMEPOINT DATA
			eFASTResults<-read.csv(paste(FILEPATH,"/",EFASTRESULTFILENAME_FULL,sep=""),header=T)
			
			TIMERESULT<-data.frame(hour,t(eFASTResults[,MEASURELABEL]))
			SiMEASURERESULTS<-rbind(SiMEASURERESULTS,TIMERESULT)				
		}
		
		colnames(SiMEASURERESULTS)<-c(TIMEPOINTSCALE,PARAMETERS)
		
		
		# PLOT THE GRAPH
		GRAPHFILE <- paste(FILEPATH,"/",MEASURE,"_OT.pdf",sep="")
		pdf(GRAPHFILE,width=7,height=7.8)
		
		GRAPHTITLE<-paste("eFAST First Order Sensitivity Indexes Over Simulation Time\nCell Response Measure: ",MEASURE,sep="")
		
		plot(TIMEPOINTS,SiMEASURERESULTS[,2],main=GRAPHTITLE,type="o",lty=1,ylim=c(0,1),pch=1,xaxt="n",xlab=TIMEPOINTSCALE,ylab="eFAST First-Order Sensitivity Index (Si)")
		# -1 TO EXCLUDE DUMMY
		for(l in 2:length(PARAMETERS)-1)
		{
			lines(TIMEPOINTS,SiMEASURERESULTS[,l+1],type="o",lty=5,pch=l)
		}
		
		axis(1,at=seq(as.numeric(min(TIMEPOINTS)),as.numeric(max(TIMEPOINTS)),by=as.numeric(max(TIMEPOINTS))/length(TIMEPOINTS)))
		legend("topleft", inset=.0,title="Parameter",PARAMETERS[1:length(PARAMETERS)-1], pch=1:length(PARAMETERS)-1,cex=0.75)
		dev.off()
	}
}


