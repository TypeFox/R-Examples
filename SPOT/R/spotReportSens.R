###################################################################################################
#'  Sensitivity Report Helper Function
#' 
#' Helper function for the sensitivity report
#'
#' @param B best solution column
#' @param fit random forest fit
#' @param roi internal parameter for the initial region of interest, or: \code{spotConfig$alg.roi}
#' @param nsens number of parameters estimated with sensitivity
#' @seealso  \code{\link{SPOT}} \code{\link{spot}} \code{\link{spotStepReport}} \code{\link{spotReportSens}} 
#' @export
#' @keywords internal
################################################################################################### 
spotReportSensY <- function(B,fit,roi,nsens) {	
  	Y <- NULL
  	for (k in 1:length(B)) {
		BV <- B;
		BV[,k] <- seq(roi[k,"lower"],roi[k,"upper"],length.out=nsens)
		Y <- data.frame(cbind(Y,predict(fit,BV)))
		names(Y)[length(Y)] <- names(B)[k]
    }	
    Y
}

###################################################################################################
#' Sensitivity Report 
#' 
#' Function to generate a report with sensitivity plot.
#'
#' The sensitivity curves are based on a metamodel which is a random forest with 100 trees 
#' fitted to the result points from RES-file. The plot contains: 
#'    x-axis:   ROI for each parameter normalized to [-1,1]
#'    y-axis: 
#'
#' @param spotConfig the configuration list of all spot parameters
#' @seealso  \code{\link{SPOT}} \code{\link{spot}} \code{\link{spotStepReport}} 
#' @export
################################################################################################### 
spotReportSens <- function(spotConfig) {		
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotReportSens")
	rawB <- spotGetRawDataMatrixB(spotConfig)
	spotPrint(spotConfig$io.verbosity,1,summary(rawB))
	mergedData <- spotPrepareData(spotConfig)
	#rawD<-spotGetRawResData(spotConfig);  # raw data with CONFIG, needed for sd(BestSolution)      MZ: replaced with if/else for fileMode
	if(spotConfig$spot.fileMode) {
		rawD <- spotGetRawResData(spotConfig)$rawD
	}else{
		rawD=spotConfig$alg.currentResult
	}	
	spotConfig=spotWriteBest(mergedData, spotConfig);
	C1=spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest),]       # WK: added ","
	xNames <- row.names(spotConfig$alg.roi)
	yNames <- spotConfig$alg.resultColumn
	B <- NULL 
	nsens=20             # number of points along the normalized ROI range
	for (i in 1:nsens) {
	   # replicate best solution (row 1), but cut away 1st column (Y):   
	   B <- rbind(B,data.frame(C1[1,xNames]));
	}
	names(B)=xNames; 	
	fit <- randomForest(rawB[xNames], rawB[[yNames]], ntree=100)		
	#fit <- rpart(y ~ ., data= rawB)
	rwb <- cbind(spotConfig$alg.roi,t(B[1,]))     # rwb: roi with 'BEST' column
	names(rwb)[length(rwb)] <- "BEST"
	Y <- spotReportSensY(B,fit,spotConfig$alg.roi,nsens)
	# scale each ROI to the normalized ROI range [-1,+1]:
	X=seq(-1,1,length.out=nsens)	
	# XP is the location of the SPOT best point for each parameter in the normalized ROI range
	XP = (rwb$BEST-rwb$lower)/(rwb$upper-rwb$lower)*2-1
	XP=rbind(XP,XP);   # matpoints needs at least two rows to plot each column in a different color
	YP = min(Y);  YP=rbind(YP,YP)
#palette("default")
	spotWriteLines(spotConfig$io.verbosity,1," ")
	spotWriteLines(spotConfig$io.verbosity,1,"Sensitivity plot for this ROI:")
	spotPrint(spotConfig$io.verbosity,1,rwb)
	spotWriteLines(spotConfig$io.verbosity,1," ")
	spotWriteLines(spotConfig$io.verbosity,1,paste("Best solution found with ",nrow(rawB)," evaluations:",sep=""))
	spotPrint(spotConfig$io.verbosity,1,C1[1,])
	spotWriteLines(spotConfig$io.verbosity,1," ")
	spotWriteLines(spotConfig$io.verbosity,1,"Standard deviation of best solution:")
	y<-rawD[rawD$CONFIG==C1[1,"CONFIG"],spotConfig$alg.resultColumn]
	spotWriteLines(spotConfig$io.verbosity,1,paste(mean(y),"+-",sd(y)))
	#finally the plot commands, both for screen and pdf file
  	if(spotConfig$report.io.pdf==TRUE){ #if pdf should be created
		pdf(spotConfig$io.pdfFileName) #start pdf creation
		matplot(X,Y,type="l",lwd=rep(3,ncol(Y)),cex.axis=1.5,cex.lab=1.5,col=1:ncol(Y),xlab="normalized ROI",main=spotConfig$userConfFileName) 
		matpoints(XP,YP,pch=rep(21,ncol(Y)),bg=1:ncol(Y),cex=2)	
		legend("topleft",legend=names(Y),lwd=rep(2,ncol(Y)),lty=1:ncol(Y),col=1:ncol(Y),text.col=1:ncol(Y))	
		dev.off() #close pdf device
	}
	if(spotConfig$report.io.screen==TRUE && spotConfig$io.verbosity>0) #if graphic should be on screen
	{
		dev.new()
		matplot(X,Y,type="l",lwd=rep(3,ncol(Y)),cex.axis=1.5,cex.lab=1.5,col=1:ncol(Y),xlab="normalized ROI",main=spotConfig$userConfFileName) 
		matpoints(XP,YP,pch=rep(21,ncol(Y)),bg=1:ncol(Y),cex=2)	
		legend("topleft",legend=names(Y),lwd=rep(2,ncol(Y)),lty=1:ncol(Y),col=1:ncol(Y),text.col=1:ncol(Y))
	}
	spotWriteLines(spotConfig$io.verbosity,2,"\n Leaving spotReportSens")
	spotConfig
}


