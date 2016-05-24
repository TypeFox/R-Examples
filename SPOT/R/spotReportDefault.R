###################################################################################################
#' Default Report
#' 
#' Function generates a simple report.
#'
#' This function is used when no report function is requested by the user. It creates some text output
#' and also draws a tree, printing it to screen or pdf. If the \code{report.io.pdf} setting is TRUE
#' the graphic is printed to a pdf file (usually named like your .conf file, and placed in the same folder)
#' if \code{report.io.screen} is set TRUE the graphic is printed to the screen. Both can be FALSE or TRUE
#' at the same time. If the user does not specify those values, the defaults will be used as shown in
#' \code{\link{spotGetOptions}}, which means there will be only screen output, and no pdf.
#'
#' In case of multi objective optimization this function will report the hyper volume indicator, and write the Pareto front and the Pareto set to
#' spotConfig$mco.val and spotConfig$mco.par. spotReportDefault is currently the only recommendable report function for multi objective SPOT, besides
#' custom functions created by users.
#' 
#' @param spotConfig the configuration list of all spot parameters
#' @return list spotConfig with changed values
#' @seealso  \code{\link{spot}}, \code{\link{spotStepReport}} 
#' @export
###################################################################################################
spotReportDefault <- function(spotConfig) {	
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotReportDefault")
	rawB <- spotGetRawDataMatrixB(spotConfig)
	spotPrint(spotConfig$io.verbosity,1,summary(rawB))
	mergedData <- spotPrepareData(spotConfig)	
	if(length(spotConfig$alg.resultColumn)==1){	
		spotConfig=spotWriteBest(mergedData, spotConfig)
		C1=spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest),]
		#cat(sprintf("\n Best solution found with %d evaluations:\n",nrow(rawB)));
		spotWriteLines(spotConfig$io.verbosity,1," ")
		spotPrint(spotConfig$io.verbosity,1,paste("Best solution found with ",nrow(rawB)," evaluations:",sep=""))
		spotPrint(spotConfig$io.verbosity,1,C1)		
		fit.tree <- rpart( paste(spotConfig$alg.resultColumn,"~",row.names(spotConfig$alg.roi)), data= rawB)
		if (!is.null(fit.tree$splits)){
			if(spotConfig$report.io.pdf==TRUE){ #if pdf should be created
				pdf(spotConfig$io.pdfFileName) #start pdf creation
				spotPlotBst(spotConfig)
				par(mfrow=c(1,1))				
				plot(fit.tree)
				text(fit.tree)
				dev.off() #close pdf device
			}
			if(spotConfig$report.io.screen==TRUE && spotConfig$io.verbosity>0) #if graphic should be on screen
			{
				dev.new()
				par(mfrow=c(1,1), xpd=NA)
				plot(fit.tree)
				text(fit.tree)
			}
		}
	}else{#else do multi criteria report 
		allPoints<-t(mergedData$mergedY)
		frontPoints<-nondominated_points(allPoints)
		if(is.null(spotConfig$mco.refPoint)){
			vol<-dominated_hypervolume(frontPoints)
		}else{
			vol<-dominated_hypervolume(frontPoints,ref=spotConfig$mco.refPoint)
		}	
		spotPrint(spotConfig$io.verbosity,1,paste("Hypervolume indicator of pareto optimal solutions: ",vol,sep=""))
		spotConfig$mco.hvolume = vol # hyper volume
		spotConfig$mco.val= t(frontPoints) #pareto front
		spotConfig$mco.par=  mergedData$x[!is_dominated(allPoints),] #pareto set
	}	
	spotConfig
}
