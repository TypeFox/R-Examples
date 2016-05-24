efast_graph_Results <-
function(RESULTS_FILE_PATH,PARAMETERS,Si,STi,errorsSi,errorsSTi,MEASURES,TIMEPOINT,TIMEPOINTSCALE)
{
	if(requireNamespace("gplots",quietly=TRUE))
	{
		colors<-c("black","grey50")

		for(MEASURE in seq(length(MEASURES)))
		{
			if(is.null(TIMEPOINT))
			{
				GRAPHFILE<-paste(RESULTS_FILE_PATH,"/",MEASURES[MEASURE],".pdf",sep="")
				GRAPHTITLE<-paste("Partitioning of Variance in Simulation Results using eFAST \n Measure: ",MEASURES[MEASURE],sep="")
			}
			else
			{
				GRAPHFILE<-paste(RESULTS_FILE_PATH,"/",MEASURES[MEASURE],"_",TIMEPOINT,".pdf",sep="")
				GRAPHTITLE<-paste("Partitioning of Variance in Simulation Results using eFAST \n Measure: ",MEASURES[MEASURE],". Timepoint: ",TIMEPOINT," ",TIMEPOINTSCALE,sep="")
			}

			pdf(GRAPHFILE)
			labelspacing<-seq(2,(length(PARAMETERS)*3),3)
		
			# DATA TO GRAPH RETRIEVES THE PARAMETERS, Si AND STi TO BE GRAPHED FROM THE MAIN RESULT SET
			dataToGraph <- data.frame(cbind(Si[,,MEASURE],STi[,,MEASURE]),check.names=FALSE)

			# CONSTRUCT THE ERROR BAR
			highSi<-dataToGraph[,1] + errorsSi[,MEASURE]
			highSTi<-dataToGraph[,2] + errorsSTi[,MEASURE]
			# COMBINE
			errorsHigh <- cbind(highSi,highSTi)

			colnames(dataToGraph)<-c("Si","STi")
			par(mar=c(9,4,4,2)+0.1)
			gplots::barplot2(t(dataToGraph),names.arg=PARAMETERS,beside=TRUE,main=GRAPHTITLE, ylim=c(0,1.0),ylab="eFAST Sensitivity",col=colors,xaxt="n",plot.ci=TRUE,ci.u=t(errorsHigh),ci.l=t(dataToGraph))
			# TEXT SIZE CONTROLLED BY CEX.AXIS
			axis(1,at=labelspacing,labels=PARAMETERS,las=2,cex.axis=0.6)
			legend("topleft",title=NULL,c("Si","STi"),fill=colors)

			dev.off()
		}
		print(paste("Graphs Output to ",RESULTS_FILE_PATH,sep=""))
	}
	else
	{
		print("Use of efast_graph_Results requires the gplots library")
	}

}

